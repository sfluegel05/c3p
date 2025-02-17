"""
Classifies: CHEBI:15889 sterol
"""
#!/usr/bin/env python
"""
Classifies: Sterol (a 3‑hydroxy steroid closely related to cholestan-3-ol)
Definition: Any molecule with a tetracyclic (or nearly so) fused ring system 
          that is cholestane‐like (≥4 rings, at least one five‐membered ring, ≥16 carbons
          and ≥75% carbon atoms in the nucleus) and with an –OH (free or glycosylated)
          attached to an atom in that nucleus.
          
Note:
  The algorithm parses the molecule, extracts all rings and their connectivity,
  chooses the largest joined ring system (candidate nucleus), checks for at least 4 rings
  and that one ring is 5-membered, then verifies that nucleus is mostly carbon.
  Finally the code looks for an –OH group (free or when linked to a sugar ring) attached
  to a nucleus carbon.
  
This is a heuristic and may fail on edge cases.
"""
from rdkit import Chem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol.
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        (bool, str): True with explanation if classified as sterol; False with reason otherwise.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get the ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # each ring is a tuple of atomic indices
    if not rings:
        return False, "No rings found in molecule"

    # Build a connectivity graph among rings (rings sharing atoms are fused).
    num_rings = len(rings)
    ring_adj = {i: set() for i in range(num_rings)}
    for i in range(num_rings):
        for j in range(i+1, num_rings):
            if set(rings[i]).intersection(rings[j]):
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

    # Choose the component having the largest number of unique atoms.
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

    # Check that the fused ring system (candidate nucleus) contains at least 4 rings.
    if len(best_component) < 4:
        return False, "Fused ring system does not contain at least 4 rings"

    # Additionally, require that at least one of the rings in the nucleus is 5-membered.
    has_5member = False
    for idx in best_component:
        if len(rings[idx]) == 5:
            has_5member = True
            break
    if not has_5member:
        return False, "Fused ring system does not contain a 5-membered ring (expected in steroids)"

    # Count carbons and total atoms in the candidate nucleus.
    carbon_count = 0
    for atom_idx in best_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:
            carbon_count += 1
    total_atoms = len(best_atoms)
    if carbon_count < 16:
        return False, f"Fused ring system has too few carbon atoms (found {carbon_count}, need >= 16)"
    if (carbon_count / total_atoms) < 0.75:
        return False, "Fused ring system contains too many heteroatoms to be a steroid nucleus"

    # FIRST: Look for a free hydroxyl group.
    # SMARTS for a free hydroxyl group: oxygen with one hydrogen.
    hydroxyl_smarts = "[OX2H]"
    hydroxyl_pattern = Chem.MolFromSmarts(hydroxyl_smarts)
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    for match in hydroxyl_matches:
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        # check that this O is attached to at least one carbon in the nucleus
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in best_atoms:
                return True, ("Molecule has a fused tetracyclic (or near‐tetracyclic) steroid nucleus "
                              "with sufficient carbon content and a free hydroxyl group attached "
                              "to the nucleus; classified as a sterol.")
    
    # SECOND: Check for a glycosylated –OH.
    # Here we look for an oxygen (not in the nucleus) attached to a nucleus carbon,
    # where that oxygen is part of a small ring having a high oxygen fraction (a sugar ring).
    for atom_idx in best_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:
            continue
        for nbr in atom.GetNeighbors():
            # We want an oxygen not included in the nucleus.
            if nbr.GetAtomicNum() != 8 or nbr.GetIdx() in best_atoms:
                continue
            # Now, check if this oxygen is part of a small ring (typically 5- or 6-membered)
            # that has a significant oxygen content.
            o_idx = nbr.GetIdx()
            for ring in rings:
                if o_idx in ring:
                    # Exclude rings that largely belong to the nucleus.
                    if len(set(ring).intersection(best_atoms)) >= len(ring):
                        continue
                    ring_size = len(ring)
                    if ring_size not in (5, 6):
                        continue
                    # Determine oxygen fraction in this candidate sugar ring.
                    oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
                    if (oxy_count / ring_size) >= 0.4:
                        return True, ("Molecule has a fused tetracyclic steroid nucleus and a glycosylated hydroxyl "
                                      "group attached to the nucleus; classified as a sterol.")
    
    return False, "No free or glycosylated hydroxyl group found attached to the fused steroid nucleus"


# Example usage:
if __name__ == "__main__":
    # Try one known sterol SMILES: (24S,25S)-cholest-5-en-3beta,24,26-triol
    test_smiles = "C[C@@H](CO)[C@@H](O)CC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"
    result, reason = is_sterol(test_smiles)
    print("Is sterol:", result)
    print("Reason:", reason)