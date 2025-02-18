"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
"""
Classifies: Phosphatidylinositol phosphate
A phosphatidylinositol phosphate (PIP, PIP2, PIP3, etc.) is a phosphoinositide that
contains a (myo-)inositol head group that is substituted by at least two phosphate groups:
one that links the inositol to the glycerol backbone and at least one extra phosphate
on the ring to distinguish it from plain phosphatidylinositol (PI). Additionally,
the molecule should contain at least two acyl chains (via ester bonds) and have a high
molecular weight.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate (PIP/PIP2/PIP3) based on its SMILES string.
    
    The classifier checks that:
      0. The molecule does not contain any negative formal charges.
      1. It contains a myo-inositol head group (based on the canonical inositol pattern).
      2. At least two distinct phosphorus atoms are directly attached to inositol (one being the linking phosphate
         which connects the inositol to a glycerol backbone, and at least one extra phosphate on the ring).
      3. At least one of the phosphate groups attached to the inositol shows “linking” behavior, meaning it is bound
         via an oxygen to a carbon atom that is not part of the inositol.
      4. The molecule contains at least two acyl ester groups (fatty acid chains attached via ester bonds) that are not
         part of any phosphate group.
      5. The molecular weight is above a lower threshold typical for these lipids.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a phosphatidylinositol phosphate, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 0. Reject molecules with any negative formal charges (to avoid deprotonated/salt forms).
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() < 0:
            return False, "Molecule has negative formal charges; likely drawn as a salt form"
    
    # 1. Locate the myo-inositol head group.
    # Use the canonical myo-inositol structure: OC1C(O)C(O)C(O)C(O)C1O
    inositol = Chem.MolFromSmiles("OC1C(O)C(O)C(O)C(O)C1O")
    inositol_matches = mol.GetSubstructMatches(inositol)
    if not inositol_matches:
        return False, "Inositol head group not found"
    # For our purposes, use the first matching fragment.
    inositol_atom_indices = set(inositol_matches[0])
    
    # 2. Collect all P atoms directly attached to any inositol oxygen.
    # We also record from which oxygen they are attached.
    p_attach = {}  # p_idx -> list of inositol oxygen indices that attach to it
    for idx in inositol_atom_indices:
        atom = mol.GetAtomWithIdx(idx)
        # Only consider oxygen atoms (the –OH groups on inositol)
        if atom.GetSymbol() != "O":
            continue
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() == "P":
                p_idx = nbr.GetIdx()
                p_attach.setdefault(p_idx, []).append(idx)
    
    # For a PIP family member we expect at least two distinct phosphate groups on inositol:
    if len(p_attach) < 2:
        return False, "No extra phosphate substituent on the inositol head found (only %d phosphate group(s) attached)" % len(p_attach)
    
    # 3. For at least one of these phosphorus atoms, check for linking phosphate behavior.
    linking_found = False
    for p_idx, o_list in p_attach.items():
        p_atom = mol.GetAtomWithIdx(p_idx)
        # Look at all oxygen neighbors of this phosphorus.
        for o_neigh in p_atom.GetNeighbors():
            # Skip if this oxygen is one of the ones attached to the inositol.
            if o_neigh.GetIdx() in inositol_atom_indices:
                continue
            if o_neigh.GetSymbol() != "O":
                continue
            # Check if this oxygen is further bound to a carbon that is not part of the inositol.
            for second_neigh in o_neigh.GetNeighbors():
                if second_neigh.GetSymbol() == "C" and second_neigh.GetIdx() not in inositol_atom_indices:
                    linking_found = True
                    break
            if linking_found:
                break
        if linking_found:
            break

    if not linking_found:
        return False, "Linking phosphate connecting inositol to glycerol backbone not found"
    
    # 4. Count acyl ester groups.
    # We define acyl ester groups as an oxygen (not involved in any phosphate) bridging a carbonyl group.
    acyl_count = 0
    seen_ester_bonds = set()
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Identify candidate pair: one atom must be oxygen and the other carbon.
        if a1.GetSymbol() == "O" and a2.GetSymbol() == "C":
            oxygen = a1
            carbon = a2
        elif a2.GetSymbol() == "O" and a1.GetSymbol() == "C":
            oxygen = a2
            carbon = a1
        else:
            continue
        
        # Exclude oxygen if it is bonded to any phosphorus (likely part of a phosphate group).
        if any(nb.GetSymbol() == "P" for nb in oxygen.GetNeighbors()):
            continue
        
        # Confirm that the carbon is part of a carbonyl (has a double bond to an oxygen aside from the ester oxygen).
        carbonyl_found = False
        for nbr in carbon.GetNeighbors():
            if nbr.GetIdx() == oxygen.GetIdx():
                continue
            if nbr.GetSymbol() != "O":
                continue
            bond_co = mol.GetBondBetweenAtoms(carbon.GetIdx(), nbr.GetIdx())
            if bond_co is not None and bond_co.GetBondType() == Chem.BondType.DOUBLE:
                carbonyl_found = True
                break
        if carbonyl_found:
            bond_id = tuple(sorted([oxygen.GetIdx(), carbon.GetIdx()]))
            if bond_id not in seen_ester_bonds:
                seen_ester_bonds.add(bond_id)
                acyl_count += 1

    if acyl_count < 2:
        return False, "Found only %d acyl ester group(s); need at least 2" % acyl_count

    # 5. Check molecular weight (typically these lipids have a high molecular weight).
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 500:
        return False, "Molecular weight too low (%.1f Da) for a phosphatidylinositol phosphate" % mw

    return True, "Molecule contains a phosphorylated inositol head (with linking and extra phosphate), and at least two acyl chains"


# Example usage (for testing):
if __name__ == "__main__":
    # Test with a known example: 1,2-dioctanoyl-sn-glycero-3-phospho-(1D-myo-inositol-3,4-bisphosphate)
    test_smiles = "CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCCC"
    result, reason = is_phosphatidylinositol_phosphate(test_smiles)
    print(result, reason)