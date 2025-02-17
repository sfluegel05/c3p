"""
Classifies: CHEBI:15889 sterol
"""
#!/usr/bin/env python
"""
Classifies: Sterol
Definition: Any 3‑hydroxy steroid whose skeleton is closely related to cholestan-3-ol.
A sterol is defined here as a molecule that:
  – Contains a fused ring system roughly corresponding to a cyclopentanoperhydrophenanthrene (steroid) nucleus,
    i.e. at least 4 fused rings with one 5-membered ring, with a high carbon fraction.
  – Contains an –OH group (either free or as part of a sugar moiety) that is directly attached to a nucleus carbon.
  
This heuristic uses a SMARTS pattern for a steroid core and then inspects the fused ring system.
Note: This is not perfect and may fail in edge cases.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol.
    
    The workflow is as follows:
      1. Parse the input SMILES.
      2. Require a positive match for a generic steroid nucleus SMARTS pattern.
      3. Extract the molecule’s ring system and choose the largest connected (fused) set.
         Check that it comprises at least 4 rings, contains at least one 5-membered ring,
         has enough carbon atoms (≥16) and a high carbon fraction (≥75%).
      4. Look for an –OH group (free or part of a small sugar ring) that is attached to a carbon in this nucleus.
    
    Args:
      smiles (str): SMILES string for the molecule.
      
    Returns:
      (bool, str): Classification decision along with a reason.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # STEP 1: Check for a steroid nucleus using a generic SMARTS pattern.
    # The pattern "C1CC2CCC3C1CCC2C3" roughly describes the cyclopentanoperhydrophenanthrene core.
    steroid_core = Chem.MolFromSmarts("C1CC2CCC3C1CCC2C3")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "Molecule does not contain a steroid nucleus (no cyclopentanoperhydrophenanthrene core detected)"
    
    # STEP 2: Analyze the ring system in the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # list of tuples, each tuple is a ring's atom indices
    if not rings:
        return False, "No rings found in molecule"
    
    # Build an adjacency graph for rings sharing atoms (i.e. fused rings).
    num_rings = len(rings)
    ring_adj = {i: set() for i in range(num_rings)}
    for i in range(num_rings):
        for j in range(i+1, num_rings):
            if set(rings[i]).intersection(rings[j]):
                ring_adj[i].add(j)
                ring_adj[j].add(i)
    
    # Identify connected components (fused systems) among the rings.
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
    
    # Choose the fused ring component with the most unique atoms.
    best_atoms = set()
    best_component = None
    for comp in components:
        comp_atoms = set()
        for idx in comp:
            comp_atoms.update(rings[idx])
        if len(comp_atoms) > len(best_atoms):
            best_atoms = comp_atoms
            best_component = comp
            
    if best_component is None or len(best_atoms) == 0:
        return False, "No fused ring system found"
    
    # Check that the fused ring system has at least 4 rings.
    if len(best_component) < 4:
        return False, "Fused ring system does not contain at least 4 rings (expected for a steroid nucleus)"
    
    # Check that at least one ring in the fused system is 5-membered.
    has_5member = any(len(rings[idx]) == 5 for idx in best_component)
    if not has_5member:
        return False, "Fused ring system does not contain a 5-membered ring (required in steroids)"
    
    # Verify the nucleus has enough carbon content.
    carbon_count = sum(1 for idx in best_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    total_atoms = len(best_atoms)
    if carbon_count < 16:
        return False, f"Fused ring system has too few carbons (found {carbon_count}, need at least 16)"
    if (carbon_count / total_atoms) < 0.75:
        return False, "Fused ring system contains too many heteroatoms to be a steroid nucleus"
    
    # Additionally, require that the steroid core SMARTS is (at least partially) located in the nucleus.
    core_matches = mol.GetSubstructMatches(steroid_core)
    found_in_nucleus = False
    for match in core_matches:
        # If a majority of the atoms in the match are in our nucleus set then we consider it valid.
        if len(set(match).intersection(best_atoms)) >= len(match) - 1:
            found_in_nucleus = True
            break
    if not found_in_nucleus:
        return False, "Steroid core pattern not found within the fused ring system"
    
    # STEP 3: Look for a hydroxyl group (free or glycosylated) attached to a carbon in the nucleus.
    # First, look for a free hydroxyl: oxygen with exactly one hydrogen.
    hydroxyl_smarts = "[OX2H]"
    hydroxyl_pattern = Chem.MolFromSmarts(hydroxyl_smarts)
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    for match in hydroxyl_matches:
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in best_atoms:
                return True, ("Molecule has a fused steroid nucleus (with ≥4 fused rings, one 5-membered, high carbon content) "
                              "and a free hydroxyl group attached to the nucleus; classified as a sterol.")
    
    # If not free, check for a glycosylated –OH:
    # Look for an oxygen (not in the nucleus) attached to a nucleus carbon and then part of a small ring (potential sugar).
    for atom_idx in best_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:
            continue
        for nbr in atom.GetNeighbors():
            # Consider only oxygen neighbors that are outside the nucleus.
            if nbr.GetAtomicNum() != 8 or nbr.GetIdx() in best_atoms:
                continue
            # See if this oxygen is part of a small ring (5 or 6 atoms) that is not largely in the nucleus.
            o_idx = nbr.GetIdx()
            for ring in rings:
                if o_idx in ring:
                    # If most of the ring is not in the nucleus, it is likely an appended sugar.
                    if len(set(ring).intersection(best_atoms)) >= len(ring):
                        continue
                    if len(ring) not in (5, 6):
                        continue
                    # Calculate oxygen fraction of the small ring.
                    oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
                    if (oxy_count / len(ring)) >= 0.4:
                        return True, ("Molecule has a fused steroid nucleus and a glycosylated hydroxyl group attached "
                                      "to the nucleus; classified as a sterol.")
    
    return False, "No free or glycosylated hydroxyl group found attached to the steroid nucleus"


# Example usage when running as a script:
if __name__ == "__main__":
    # Test with one known sterol SMILES: (24S,25S)-cholest-5-en-3beta,24,26-triol
    test_smiles = "C[C@@H](CO)[C@@H](O)CC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"
    result, reason = is_sterol(test_smiles)
    print("Is sterol:", result)
    print("Reason:", reason)