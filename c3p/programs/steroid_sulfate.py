"""
Classifies: CHEBI:16158 steroid sulfate
"""
#!/usr/bin/env python
"""
Classifies: steroid sulfate
Definition: A sulfuric ester obtained by the formal condensation of a hydroxy group 
of any steroid with sulfuric acid. We require that a valid steroid sulfate must have 
a sulfate ester group (anionic or neutral) that is attached (directly or via one intervening saturated atom) 
to a steroid nucleus. Our heuristic for a steroid nucleus is a fused polycyclic system composed solely 
of rings of size 5 or 6, with at least four rings (at least one 5-membered and three 6-membered rings).
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    
    Process:
    1. Parse the molecule.
    2. Identify a candidate steroid nucleus:
       - Collect only rings of size 5 or 6.
       - Group rings that are fused (share at least one atom).
       - In any fused group, require at least 4 rings overall 
         and (by counting ring sizes) at least one five‐membered ring and three six‐membered rings.
    3. Identify sulfate ester attachment oxygen(s):
       - Look at O atoms that are bonded (by a single bond) to an S atom that is bonded to three oxygens.
       - Exclude oxygens that are double‐bonded to S.
    4. For each candidate sulfate attachment oxygen, check if there is a connection 
       (via a path of length 1 or 2 bonds; if 2 bonds the intervening atom must be sp3 and nonaromatic)
       to any atom in the steroid nucleus.
    
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
    # Extra pass to merge indirectly fused groups
    merged = True
    while merged:
        merged = False
        new_groups = []
        while fused_groups:
            grp = fused_groups.pop()
            found = False
            for idx, other in enumerate(new_groups):
                if grp & other:
                    new_groups[idx] = other.union(grp)
                    found = True
                    merged = True
                    break
            if not found:
                new_groups.append(grp)
        fused_groups = new_groups

    # Next, for each fused group, count how many rings (from rings_5_or_6) are fully contained within it 
    # and count how many are 5-membered and 6-membered.
    steroid_candidate = None
    for group in fused_groups:
        contained_rings = [ring for ring in rings_5_or_6 if ring.issubset(group)]
        if len(contained_rings) >= 4:
            count_5 = sum(1 for ring in contained_rings if len(ring) == 5)
            count_6 = sum(1 for ring in contained_rings if len(ring) == 6)
            if count_5 >= 1 and count_6 >= 3:
                steroid_candidate = group
                break
    if steroid_candidate is None:
        return False, ("No fused tetracyclic (steroid) nucleus found "
                       "(requires at least 4 fused rings with 1 five‐membered and 3 six‐membered rings)")
    
    # For faster lookup, mark candidate nucleus atoms
    nucleus_atoms = steroid_candidate

    # ----- Step 2: Identify sulfate ester attachment oxygens -----
    candidate_oxy_atoms = []
    # We loop over all oxygen atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        # We want an oxygen that is single-bonded to an S atom.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 16:  # sulfur
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond.GetBondType() != rdchem.BondType.SINGLE:
                    continue
                # Check that this S is “sulfate-like”: it must be bonded to three oxygens 
                # (typically two via double bonds and one single-bond not back to this O)
                o_count = sum(1 for n in nbr.GetNeighbors() if n.GetAtomicNum() == 8)
                if o_count < 3:
                    continue
                # Now, ensure that our oxygen is not a double-bonded O to the S.
                # (We check the bond order between this O and S.)
                if bond.GetBondTypeAsDouble() >= 2.0:
                    continue
                candidate_oxy_atoms.append(atom.GetIdx())
    if not candidate_oxy_atoms:
        return False, "No sulfate ester attachment oxygen found"
    
    # ----- Step 3: Check connectivity from a candidate oxygen to the steroid nucleus -----
    # Allowed path lengths: 1 (directly attached) or 2 bonds (via one linker atom)
    found_attachment = False
    for o_idx in candidate_oxy_atoms:
        # For each atom in the nucleus, get the shortest path from this candidate oxygen.
        for nucleus_idx in nucleus_atoms:
            sp = Chem.GetShortestPath(mol, o_idx, nucleus_idx)
            if sp:
                # sp is a tuple of atom indices along the path.
                path_length = len(sp) - 1  # number of bonds
                if path_length == 1:
                    found_attachment = True
                    break
                elif path_length == 2:
                    # Check the linker atom (the one in between)
                    linker_idx = sp[1]
                    linker_atom = mol.GetAtomWithIdx(linker_idx)
                    # allow only if linker is sp3 and not aromatic
                    if (linker_atom.GetHybridization() == rdchem.HybridizationType.SP3
                        and not linker_atom.GetIsAromatic()):
                        found_attachment = True
                        break
                # If longer than 2, ignore.
        if found_attachment:
            break

    if not found_attachment:
        return False, "Sulfate group not attached (even via one saturated linker) to the steroid nucleus"
    
    return True, "Contains sulfate ester group attached (directly or via a short linker) to steroid nucleus (fused tetracyclic system)"

# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = [
        # True positives (should be classified as steroid sulfate):
        "[Na+].C[C@]12CC[C@H]3C(=CCc4cc(OS([O-])(=O)=O)ccc34)[C@@H]1CCC2=O",  # equilin sodium sulfate
        "C[C@]12CC[C@@H]3[C@H]([C@H]1CCC2=O)CCC4=C3C=CC(=C4)OS(=O)(=O)O",  # steroid sulfate ester (sulfuric acid ester)
        "S(O[C@H]1C[C@@]2([C@@]3([C@]([C@]4([C@](CC3)([C@H]([C@@H](O)[C@H]4O)C=C)C)[H])(CC[C@]2([C@@H](O)[C@H]1O)[H])[H])[H])C)(O)(=O)=O",  # Ptilosteroid A
        "S(O[C@H]1C[C@@]2([C@@](C3C(C4[C@](CC3)(C(=O)CC4)C)CC2)(CC1)C)[H])(O)(=O)=O",  # Etiocholanolone sulfate
        # False negatives (should be classified as steroid sulfate but missed before):
        "S(OCC(=O)[C@@H]1[C@@]2([C@]([C@]3([C@@]([C@@]4(C(=CC3)C[C@@H](O)CC4)C)(CC2)[H])[H])(CC1)[H])C)(O)(=O)=O",  # 21-hydroxypregnenolone monosulfate
        # Molecules that are not steroid sulfates:
        "O[S](=O)(=O)O",  # pure sulfate group
        "CCCC",         # non-cyclic aliphatic molecule
    ]
    for smi in test_smiles:
        result, reason = is_steroid_sulfate(smi)
        print("SMILES:", smi)
        print("Result:", result, "|", reason)
        print("-"*80)