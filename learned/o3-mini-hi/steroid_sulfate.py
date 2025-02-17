"""
Classifies: CHEBI:16158 steroid sulfate
"""
#!/usr/bin/env python
"""
Classifies: steroid sulfate
Definition: A sulfuric ester obtained by the formal condensation of a hydroxy group of any steroid with sulfuric acid.
The molecule must contain a sulfate ester group attached (directly or via the oxygen of a hydroxy substituent) to a steroid nucleus.
Heuristic for a steroid nucleus: a fused polycyclic system that contains at least four rings, with at least one 5-membered ring and three 6-membered rings.
"""

from rdkit import Chem

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    
    The classification checks:
    1. The presence of a sulfate ester group (either in anionic form like -OS(=O)(=O)[O-] or neutral -OS(=O)(=O)O).
    2. The presence of a steroid nucleus, which is heuristically defined as a fused system 
       containing at least 4 rings – with at least one 5-membered and at least three 6-membered rings.
    3. That at least one sulfate ester is attached via its oxygen (coming from the original hydroxy group)
       to an atom of the steroid nucleus or directly adjacent (i.e. the atom bonded to that oxygen is in the nucleus).
       
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a steroid sulfate; False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ------ 1. Identify sulfate ester group(s) --------
    # Define two patterns: one for anionic sulfate ester and one for neutral sulfate ester.
    sulfate_pattern_anionic = Chem.MolFromSmarts("O[S](=O)(=O)[O-]")
    sulfate_pattern_neutral = Chem.MolFromSmarts("O[S](=O)(=O)O")
    sulfate_matches = []
    if sulfate_pattern_anionic:
        sulfate_matches.extend(mol.GetSubstructMatches(sulfate_pattern_anionic))
    if sulfate_pattern_neutral:
        sulfate_matches.extend(mol.GetSubstructMatches(sulfate_pattern_neutral))
    if not sulfate_matches:
        return False, "No sulfate ester group found"
    
    # ------ 2. Identify steroid nucleus via fused rings ---------
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # tuple of tuples with atom indices
    if not atom_rings:
        return False, "No rings found; not a steroid nucleus"
    
    # Create list of sets of atom indices for each ring.
    rings = [set(ring) for ring in atom_rings]
    
    # Group rings that are fused (i.e. share at least one atom) into fused groups.
    fused_groups = []
    for ring in rings:
        added = False
        for group in fused_groups:
            if ring & group:  # if they share an atom
                group.update(ring)
                added = True
                break
        if not added:
            fused_groups.append(set(ring))
    # Sometimes groups can merge indirectly – do an extra merge pass.
    merged = True
    while merged:
        merged = False
        new_groups = []
        while fused_groups:
            grp = fused_groups.pop()
            found_merge = False
            for i, h in enumerate(new_groups):
                if grp & h:
                    new_groups[i] = h.union(grp)
                    found_merge = True
                    merged = True
                    break
            if not found_merge:
                new_groups.append(grp)
        fused_groups = new_groups

    # For each fused group, count how many of the original rings are completely contained in it
    # and record the sizes of these rings.
    candidate_group = None
    for group in fused_groups:
        contained_rings = []
        for ring in rings:
            if ring.issubset(group):
                contained_rings.append(ring)
        # Heuristic: at least 4 rings needed.
        if len(contained_rings) >= 4:
            # Count ring sizes.
            sizes = [len(ring) for ring in contained_rings]
            count_5 = sum(1 for s in sizes if s == 5)
            count_6 = sum(1 for s in sizes if s == 6)
            if count_5 >= 1 and count_6 >= 3:
                candidate_group = group
                break
    if candidate_group is None:
        return False, "No fused tetracyclic (steroid) nucleus found (require 1 five-membered and 3 six-membered rings)"

    # ------ 3. Check that a sulfate ester is attached to the steroid nucleus -----------
    # For the sulfate match, the pattern (O[S](=O)(=O)[O-] or O[S](=O)(=O)O) returns four atoms.
    # Atom index 0 is the oxygen that originally comes from the hydroxy group.
    attached = False
    for match in sulfate_matches:
        # match is a tuple like (o_idx, s_idx, o1_idx, o2_idx)
        att_o_idx = match[0]  # oxygen that should be coming from the steroid hydroxy group.
        att_o = mol.GetAtomWithIdx(att_o_idx)
        # Check neighbors of this oxygen that are not sulfur; one of these should be in (or adjacent to) the steroid nucleus.
        for nbr in att_o.GetNeighbors():
            # Skip the sulfate sulfur.
            if nbr.GetSymbol() == "S":
                continue
            # If the neighbor is in the candidate fused ring group, we consider it attached.
            if nbr.GetIdx() in candidate_group:
                attached = True
                break
            # Alternatively, if the neighbor is not itself in the ring but is connected to an atom in the nucleus (e.g. via a single bond), allow that.
            for nn in nbr.GetNeighbors():
                if nn.GetIdx() in candidate_group:
                    attached = True
                    break
            if attached:
                break
        if attached:
            break
    if not attached:
        return False, "Sulfate group not attached to the steroid nucleus"

    return True, "Contains sulfate ester group attached to steroid nucleus (fused tetracyclic system)"

# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = [
        "[Na+].C[C@]12CC[C@H]3C(=CCc4cc(OS([O-])(=O)=O)ccc34)[C@@H]1CCC2=O",  # equilin sodium sulfate (TP)
        "O[S](=O)(=O)O",  # just a sulfate group molecule (should fail, no steroid)
        "CCCC",         # totally non-steroidal 
        # add additional SMILES from the outcomes if desired.
    ]
    for smi in test_smiles:
        result, reason = is_steroid_sulfate(smi)
        print(smi, "->", result, "|", reason)