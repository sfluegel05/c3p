"""
Classifies: CHEBI:15889 sterol
"""
#!/usr/bin/env python
"""
Classifies a sterol as defined by:
  "Any 3‑hydroxy steroid whose skeleton is closely related to cholestan-3-ol 
   (additional carbon atoms may be present in the side chain)."

Our improved approach:
  1. Parse the SMILES string.
  2. Get all rings in the molecule and build a connectivity graph between rings.
  3. Identify the largest fused ring system (by unique atom count).
  4. Require that the fused ring system consists of at least 4 rings, includes at least 16 carbons 
     and is mostly carbon (>=75% of atoms).
  5. Look for a hydroxyl group directly attached to a nucleus carbon. If not free, check whether 
     the oxygen is linked to a sugar‐like ring (i.e. a glycosylated hydroxyl typical of sterol glycosides).
"""

from rdkit import Chem
import math

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol.
    A sterol is defined as a 3-hydroxy steroid whose fused (typically tetracyclic)
    ring system is closely related to cholestan-3-ol. In most sterols the steroid nucleus 
    is formed by three six-membered rings and one five-membered ring, and a hydroxyl group 
    is attached (directly or via a glycosidic bond) to one of the ring carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a sterol, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # each ring is a tuple of atomic indices
    if not rings:
        return False, "No rings found in molecule"
    
    # Build a connectivity graph between rings.
    num_rings = len(rings)
    ring_adj = {i: set() for i in range(num_rings)}
    for i in range(num_rings):
        for j in range(i+1, num_rings):
            if set(rings[i]).intersection(rings[j]):  # shared atoms => fused rings
                ring_adj[i].add(j)
                ring_adj[j].add(i)
    
    # Identify connected components among rings.
    visited = set()
    components = []
    for i in range(num_rings):
        if i in visited:
            continue
        comp = set()
        stack = [i]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            comp.add(node)
            stack.extend(ring_adj[node] - visited)
        components.append(comp)
    
    # Choose the component that has the largest number of unique atoms (candidate steroid nucleus).
    best_atoms = set()
    best_component = None
    for comp in components:
        comp_atoms = set()
        for idx in comp:
            comp_atoms.update(rings[idx])
        if len(comp_atoms) > len(best_atoms):
            best_atoms = comp_atoms
            best_component = comp
    
    if best_component is None:
        return False, "No fused ring system found"

    # In a typical steroid nucleus, there are 4 fused rings.
    if len(best_component) < 4:
        return False, "Fused ring system does not contain at least 4 rings"
    
    # Count carbons and total atoms in the candidate nucleus.
    carbon_count = 0
    for atom_idx in best_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:
            carbon_count += 1
    total_atoms_in_nucleus = len(best_atoms)
    
    if carbon_count < 16:
        return False, f"Fused ring system has too few carbons (found {carbon_count}, need >= 16)"
    if (carbon_count / total_atoms_in_nucleus) < 0.75:
        return False, "Fused ring system contains too many heteroatoms to be a steroid nucleus"
    
    # Function to search for free hydroxyl groups (pattern: oxygen with one hydrogen)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    # Check for a free hydroxyl group directly attached to a nucleus carbon.
    for match in hydroxyl_matches:
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        # A free -OH group should be attached to at least one carbon.
        for neighbor in o_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() in best_atoms:
                return True, ("Molecule has a fused tetracyclic steroid nucleus (>=4 rings, >=16 carbons, mostly carbon) "
                              "and a free hydroxyl group attached to the nucleus; classified as a sterol.")
    
    # If no free -OH is found, check for glycosylated hydroxyl: an oxygen attached to a nucleus carbon,
    # which in turn is bonded to a sugar-like ring (small ring with multiple oxygens).
    # Get all rings (we already have them as 'rings').
    for atom_idx in best_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        # Only consider carbon atoms; typically the hydroxyl should be linked to a C in the nucleus.
        if atom.GetAtomicNum() != 6:
            continue
        for neighbor in atom.GetNeighbors():
            # Look for an oxygen not in the nucleus.
            if neighbor.GetAtomicNum() != 8 or neighbor.GetIdx() in best_atoms:
                continue
            # Now check if that oxygen is part of a small ring with a high oxygen content (possible sugar ring).
            o_idx = neighbor.GetIdx()
            for ring in rings:
                # Consider ring if it does not largely overlap with the nucleus.
                if o_idx in ring and len(set(ring).intersection(best_atoms)) < len(ring):
                    ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
                    ring_size = len(ring_atoms)
                    # Only consider rings of size 5 or 6 as typical sugar rings.
                    if ring_size not in (5,6):
                        continue
                    oxygen_in_ring = sum(1 for a in ring_atoms if a.GetAtomicNum() == 8)
                    # If the fraction of oxygens is high enough (e.g. >= 40%), count it as sugar-like.
                    if (oxygen_in_ring/ring_size) >= 0.4:
                        return True, ("Molecule has a fused tetracyclic steroid nucleus and a glycosylated hydroxyl "
                                      "attached via an oxygen to the nucleus; classified as a sterol.")
    
    return False, "No hydroxyl (free or glycosylated) found attached to the fused steroid nucleus"

# Example usage:
if __name__ == "__main__":
    # Test using one of the example SMILES (a known sterol)
    test_smiles = "C[C@@H](CO)[C@@H](O)CC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"
    result, reason = is_sterol(test_smiles)
    print("Is sterol:", result)
    print("Reason:", reason)