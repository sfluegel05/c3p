"""
Classifies: CHEBI:16158 steroid sulfate
"""
#!/usr/bin/env python
"""
Classifies: steroid sulfate
Definition: A sulfuric ester obtained by the formal condensation of a hydroxy group 
of any steroid with sulfuric acid. For our purposes a valid steroid sulfate must have a sulfate ester group (either anionic or neutral) 
that is attached either directly or via one intervening saturated (sp3 non‐aromatic) linker to a steroid nucleus. 
Our heuristic for a steroid nucleus is a fused polycyclic system composed solely
of 5‐ and 6‐membered rings with at least four rings (and with at least one 5‐membered ring and three 6‐membered rings).
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    
    Process:
    1. Parse the molecule.
    2. Identify a candidate steroid nucleus:
       - Extract all rings of size 5 or 6.
       - Group rings that are fused (share at least one atom, merging also indirect fusion).
       - Accept if any fused group contains at least 4 rings overall with at least 1 five‐membered ring and 3 six‐membered rings.
    3. Identify sulfate ester attachment oxygens:
       - Look for O atoms that are single-bonded to a S atom.
       - The sulfur must be “sulfate‐like” (bonded to at least three oxygens).
    4. For each candidate attachment oxygen, check connectivity:
       - If the oxygen’s neighbor (other than the sulfate S) is in the nucleus, accept.
       - Otherwise, if that neighbor is sp3, non‐aromatic and has any neighbor that belongs to the steroid nucleus, accept.
       
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a steroid sulfate; False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ----- Step 1: Identify candidate steroid nucleus -----
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    if not all_rings:
        return False, "No rings found; not a steroid nucleus"
    
    # Consider only rings of size 5 or 6.
    rings_5_or_6 = [set(ring) for ring in all_rings if len(ring) in (5, 6)]
    if not rings_5_or_6:
        return False, "No 5- or 6-membered rings found; not a steroid nucleus"
    
    # Group fused rings (rings sharing at least one atom)
    fused_groups = []
    for ring in rings_5_or_6:
        merged = False
        for group in fused_groups:
            if ring & group:  # if they share at least one atom
                group.update(ring)
                merged = True
                break
        if not merged:
            fused_groups.append(set(ring))
    # Merge indirectly fused groups in a second pass.
    merged_any = True
    while merged_any:
        merged_any = False
        new_groups = []
        while fused_groups:
            grp = fused_groups.pop()
            found = False
            for idx, other in enumerate(new_groups):
                if grp & other:
                    new_groups[idx] = other.union(grp)
                    found = True
                    merged_any = True
                    break
            if not found:
                new_groups.append(grp)
        fused_groups = new_groups

    # Next, from all rings, count which fused groups yield a steroid nucleus prediction.
    steroid_candidate = None
    for group in fused_groups:
        # Count rings (from our rings_5_or_6) that are fully contained in this group.
        contained_rings = [ring for ring in rings_5_or_6 if ring.issubset(group)]
        if len(contained_rings) >= 4:
            count_5 = sum(1 for ring in contained_rings if len(ring)==5)
            count_6 = sum(1 for ring in contained_rings if len(ring)==6)
            if count_5 >= 1 and count_6 >= 3:
                steroid_candidate = group
                break
    if steroid_candidate is None:
        return False, ("No fused tetracyclic (steroid) nucleus found "
                       "(requires at least 4 fused rings with at least 1 five-membered and 3 six-membered rings)")
    nucleus_atoms = steroid_candidate  # set of atom indices in the nucleus
    
    # ----- Step 2: Identify sulfate ester attachment oxygens -----
    candidate_oxy_indices = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue  # not oxygen
        # Check if bonded to a sulfur
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 16:  # sulfur
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                # must be a single bond
                if bond.GetBondType() != rdchem.BondType.SINGLE:
                    continue
                # Check that this sulfur is “sulfate-like”
                # Count oxygen neighbors of the sulfur.
                oxy_count = 0
                for snbr in nbr.GetNeighbors():
                    if snbr.GetAtomicNum() == 8:
                        oxy_count += 1
                if oxy_count < 3:
                    continue
                # Make sure our oxygen is not one that is double-bonded to S.
                if bond.GetBondTypeAsDouble() >= 2.0:
                    continue
                candidate_oxy_indices.append(atom.GetIdx())
                break  # no need to check additional neighbors for this oxygen

    if not candidate_oxy_indices:
        return False, "No sulfate ester attachment oxygen found"
    
    # ----- Step 3: Check connectivity from candidate oxygen to the steroid nucleus -----
    attached_found = False
    for oxy_idx in candidate_oxy_indices:
        oxy_atom = mol.GetAtomWithIdx(oxy_idx)
        # Look at neighbors other than the sulfur (usually there is only one such neighbor)
        for nbr in oxy_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 16:
                continue  # skip the sulfate sulfur neighbor
            # Case 1: Direct attachment of the oxygen to an atom in the nucleus.
            if nbr.GetIdx() in nucleus_atoms:
                attached_found = True
                break
            # Case 2: One-linker attachment.
            # The oxygen is attached to a linker atom that is not in the nucleus but that linker must be sp3 and non-aromatic.
            # Then that linker should be attached to at least one nucleus atom.
            if nbr.GetHybridization() != rdchem.HybridizationType.SP3 or nbr.GetIsAromatic():
                continue
            for nn in nbr.GetNeighbors():
                if nn.GetIdx() == oxy_idx:
                    continue
                if nn.GetIdx() in nucleus_atoms:
                    attached_found = True
                    break
            if attached_found:
                break
        if attached_found:
            break

    if not attached_found:
        return False, "Sulfate group not attached (even via one saturated linker) to the steroid nucleus"
    
    return True, "Contains sulfate ester group attached (directly or via one saturated linker) to steroid nucleus (fused tetracyclic system)"

# Example usage (testing with sample SMILES):
if __name__ == "__main__":
    test_smiles = [
        # True positives:
        "[Na+].C[C@]12CC[C@H]3C(=CCc4cc(OS([O-])(=O)=O)ccc34)[C@@H]1CCC2=O",  # equilin sodium sulfate
        "C[C@]12CC[C@@H]3[C@H]([C@H]1CCC2=O)CCC4=C3C=CC(=C4)OS(=O)(=O)O",  # sulfuric acid ester example
        "S(O[C@H]1C[C@@]2([C@@]3([C@]([C@]4([C@](CC3)([C@H]([C@@H](O)[C@H]4O)C=C)C)[H])(CC[C@]2([C@@H](O)[C@H]1O)[H])[H])[H])C)(O)(=O)=O",  # Ptilosteroid A
        # False negatives (should be steroid sulfate):
        "S(OCC(=O)[C@@H]1[C@@]2([C@]([C@]3([C@@]([C@@]4(C(=CC3)C[C@@H](O)CC4)C)(CC2)[H])[H])(CC1)[H])C)(O)(=O)=O",  # 21-hydroxypregnenolone monosulfate
        # Molecules that are NOT steroid sulfates (should be rejected):
        "O[S](=O)(=O)O",  # pure sulfate group
        "CCCC",         # non-cyclic aliphatic
    ]
    for smi in test_smiles:
        result, reason = is_steroid_sulfate(smi)
        print("SMILES:", smi)
        print("Result:", result, "|", reason)
        print("-"*80)