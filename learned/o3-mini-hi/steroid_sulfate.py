"""
Classifies: CHEBI:16158 steroid sulfate
"""
#!/usr/bin/env python
"""
Classifies: steroid sulfate
Definition: A sulfuric ester obtained by the formal condensation of a hydroxy group of any steroid with sulfuric acid.
"""

from rdkit import Chem

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    
    The classification is based on two criteria:
    1. The presence of a sulfate ester group (–OS(=O)(=O)O, or its ionic form).
    2. The presence of a steroid nucleus, which we heuristically define as a set of
       at least four fused rings (usually three six-membered rings and one five-membered ring)
       that are interconnected. Additionally, one of the sulfate groups should be attached
       to one of the atoms within that fused ring system (i.e. the steroid framework).
       
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a steroid sulfate; False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the sulfate ester group SMARTS.
    # This pattern matches the O atom that is bound to sulfur which in turn is bound to two =O groups.
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)")
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)
    if not sulfate_matches:
        return False, "No sulfate ester group found"
    
    # Identify all rings in the molecule.
    ring_info = mol.GetRingInfo()
    rings = [set(ring) for ring in ring_info.AtomRings()]
    if not rings:
        return False, "No rings found in the molecule, not a steroid sulfate"
    
    # Group rings that are fused (i.e. share at least one atom) into connected components.
    groups = []
    for ring in rings:
        added = False
        for grp in groups:
            # If this ring shares at least one atom with the group, merge it.
            if ring & grp:
                grp.update(ring)
                added = True
                break
        if not added:
            groups.append(set(ring))
    
    # Sometimes groups may merge further if rings are connected indirectly.
    merged = True
    while merged:
        merged = False
        new_groups = []
        while groups:
            g = groups.pop()
            # Try to merge with any group already in new_groups.
            found_merge = False
            for i, h in enumerate(new_groups):
                if g & h:
                    new_groups[i] = h.union(g)
                    found_merge = True
                    merged = True
                    break
            if not found_merge:
                new_groups.append(g)
        groups = new_groups
    
    # Heuristic: Look for a fused ring group that contains at least four rings.
    # We again use the ring information: count how many of the original rings are fully contained in each group.
    candidate = None
    for grp in groups:
        count = 0
        for ring in rings:
            if ring.issubset(grp):
                count += 1
        if count >= 4:
            candidate = grp
            break
    if candidate is None:
        return False, "No fused tetracyclic (steroid) nucleus found"
    
    # Now, check that at least one sulfate group is directly attached to the steroid nucleus.
    # For each sulfate match, the pattern returns the index of the O atom in "OS(=O)(=O)".
    attached = False
    for match in sulfate_matches:
        # In our sulfate_pattern, atom index 0 corresponds to the oxygen that is bound to sulfur.
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        # Get neighbors (the sulfate oxygen is bonded to sulfur and to one more atom).
        for nbr in o_atom.GetNeighbors():
            # Skip if the neighbor is sulfur (which is part of the sulfate fragment).
            if nbr.GetSymbol() == "S":
                continue
            # If the neighboring atom is part of the fused ring candidate, we assume attachment.
            if nbr.GetIdx() in candidate:
                attached = True
                break
        if attached:
            break
    if not attached:
        return False, "Sulfate group not attached to the steroid nucleus"

    return True, "Contains a sulfate ester group attached to a steroid nucleus (fused tetracyclic system)"

# Example usage and tests (uncomment for testing):
# test_smiles = [
#     "[Na+].C[C@]12CC[C@H]3C(=CCc4cc(OS([O-])(=O)=O)ccc34)[C@@H]1CCC2=O",  # equilin sodium sulfate
#     "C[C@]12CC[C@@H]3[C@H]([C@H]1CCC2=O)CCC4=C3C=CC(=C4)OS(=O)(=O)O",  # another example
#     "CCCC",  # non-steroid
# ]
# for smi in test_smiles:
#     result, reason = is_steroid_sulfate(smi)
#     print(smi, "->", result, reason)