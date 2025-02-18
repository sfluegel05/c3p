"""
Classifies: CHEBI:16158 steroid sulfate
"""
#!/usr/bin/env python
"""
Classifies: steroid sulfate
Definition: A sulfuric ester obtained by the formal condensation of a hydroxy group of any steroid with sulfuric acid.
A valid steroid sulfate must have a sulfate ester group (in anionic or neutral form)
that is attached (directly or via one intervening saturated atom) to a steroid nucleus.
Heuristic for a steroid nucleus: a fused polycyclic system with at least 4 rings,
including at least one 5-membered and three 6-membered rings.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    
    Steps:
    1. Look for sulfate ester group (either anionic or neutral).
    2. Identify a steroid nucleus via a fused group of rings:
       - at least four rings,
       - at least one 5-membered and three 6-membered rings.
    3. Verify that one sulfate ester’s "attachment oxygen" (i.e. the oxygen from the -OH that underwent condensation)
       is connected (directly or via one intervening sp3 saturated atom) to an atom that is part of the steroid nucleus.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a steroid sulfate; False otherwise.
        str: Explanation of classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ---- Step 1: Identify sulfate ester group(s) ----
    # We define two SMARTS: one for anionic sulfate ester, one for neutral sulfate ester.
    sulfate_smarts_a = "O[S](=O)(=O)[O-]"  # anionic form
    sulfate_smarts_n = "O[S](=O)(=O)O"     # neutral form
    sulfate_pattern_a = Chem.MolFromSmarts(sulfate_smarts_a)
    sulfate_pattern_n = Chem.MolFromSmarts(sulfate_smarts_n)
    sulfate_matches = []
    if sulfate_pattern_a:
        sulfate_matches.extend(mol.GetSubstructMatches(sulfate_pattern_a))
    if sulfate_pattern_n:
        sulfate_matches.extend(mol.GetSubstructMatches(sulfate_pattern_n))
    if not sulfate_matches:
        return False, "No sulfate ester group found"
    
    # ---- Step 2: Identify steroid nucleus via fused rings ----
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # tuple of tuples of atom indices
    if not atom_rings:
        return False, "No rings found; not a steroid nucleus"
    
    # Convert each ring into a set of atom indices.
    rings = [set(ring) for ring in atom_rings]
    
    # Group rings that share atoms (fused) into fused groups.
    fused_groups = []
    for ring in rings:
        added = False
        for group in fused_groups:
            if ring & group:   # any common atom means fused
                group.update(ring)
                added = True
                break
        if not added:
            fused_groups.append(set(ring))
    # Do an extra pass for indirect merging.
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

    # For each fused group, count how many rings (from our original ring list) are completely contained in it.
    candidate_group = None
    for group in fused_groups:
        contained_rings = [ring for ring in rings if ring.issubset(group)]
        if len(contained_rings) >= 4:  # at least four rings present
            # Check ring sizes
            sizes = [len(ring) for ring in contained_rings]
            count_5 = sum(1 for s in sizes if s == 5)
            count_6 = sum(1 for s in sizes if s == 6)
            if count_5 >= 1 and count_6 >= 3:
                candidate_group = group
                break
    if candidate_group is None:
        return False, "No fused tetracyclic (steroid) nucleus found (requires 1 five-membered and 3 six-membered rings)"
    
    # ---- Step 3: Check that a sulfate ester is attached to the steroid nucleus ----
    # For each sulfate match, we assume that the first atom in the match (index 0) is the oxygen originating from the hydroxy group.
    # We check if this oxygen is (a) sp3 (typical for a former -OH) and (b) is connected (directly or within two bonds) to an atom in candidate_group.
    attached = False
    for match in sulfate_matches:
        # match example: (o_idx, s_idx, o1_idx, o2_idx)
        att_o_idx = match[0]
        att_o = mol.GetAtomWithIdx(att_o_idx)
        # Optional: check that the oxygen is not in an aromatic system (suggesting it is an -OH fragment)
        if att_o.GetIsAromatic():
            continue  # skip if aromatic
        
        # First check: direct neighbors
        for nbr in att_o.GetNeighbors():
            if nbr.GetSymbol() == "S":
                continue  # ignore the sulfur in the sulfate
            if nbr.GetIdx() in candidate_group:
                attached = True
                break
            # Second check: neighbor-of-neighbor (only if the connecting neighbor is sp3 saturated, indicating a linker)
            if nbr.GetHybridization() == rdchem.HybridizationType.SP3:
                for nn in nbr.GetNeighbors():
                    if nn.GetIdx() == att_o_idx:
                        continue
                    if nn.GetIdx() in candidate_group:
                        attached = True
                        break
            if attached:
                break
        if attached:
            break
    if not attached:
        return False, "Sulfate group not attached (even via a short saturated linker) to the steroid nucleus"
    
    return True, "Contains sulfate ester group attached (directly or via a short linker) to steroid nucleus (fused tetracyclic system)"

# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = [
        # True positives (should be classified as steroid sulfate):
        "[Na+].C[C@]12CC[C@H]3C(=CCc4cc(OS([O-])(=O)=O)ccc34)[C@@H]1CCC2=O",  # equilin sodium sulfate
        "C[C@]12CC[C@@H]3[C@H]([C@H]1CCC2=O)CCC4=C3C=CC(=C4)OS(=O)(=O)O",  # sulfuric acid ester of steroid
        "S(O[C@H]1C[C@@]2([C@@]3([C@]([C@]4([C@](CC3)([C@H]([C@@H](O)[C@H]4O)C=C)C)[H])(CC[C@]2([C@@H](O)[C@H]1O)[H])[H])[H])C)(O)(=O)=O",  # Ptilosteroid A
        # False negatives (should be classified as steroid sulfate but were missed before):
        "S(OCC(=O)[C@@H]1[C@@]2([C@]([C@]3([C@@]([C@@]4(C(=CC3)C[C@@H](O)CC4)C)(CC2)[H])[H])(CC1)[H])C)(O)(=O)=O",  # 21-hydroxypregnenolone monosulfate
        # Non-steroidal examples:
        "O[S](=O)(=O)O",  # pure sulfate group
        "CCCC",         # non-cyclic aliphatic
    ]
    for smi in test_smiles:
        result, reason = is_steroid_sulfate(smi)
        print("SMILES:", smi)
        print("Result:", result, "|", reason)
        print("-"*80)