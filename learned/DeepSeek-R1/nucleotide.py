"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: CHEBI:36976 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    A nucleotide consists of a nucleoside (sugar + base) with at least one phosphate group
    attached to the 3' or 5' hydroxyl group of the sugar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Improved nucleoside detection: sugar + nitrogenous base
    # 1. Find furanose sugar pattern (C1OCC(O)CO with stereochemistry)
    sugar_pattern = Chem.MolFromSmarts("[C@H]1([C@H](O)[C@@H](O)C(O)[CH]1)")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No ribose-like sugar detected"

    # 2. Check for nitrogenous base attached to C1' (first carbon in sugar)
    base_attached = False
    for sugar_match in sugar_matches:
        sugar_atoms = list(sugar_match)
        c1 = sugar_atoms[0]  # First atom in sugar pattern is C1'
        for neighbor in mol.GetAtomWithIdx(c1).GetNeighbors():
            if neighbor.GetAtomicNum() == 7 and neighbor.GetIsAromatic():
                base_attached = True
                break
        if base_attached:
            break
    if not base_attached:
        return False, "No nitrogenous base attached to sugar"

    # 3. Check phosphate groups attached to 3' or 5' positions
    # 5' is first carbon in pattern (C1'), 3' is third carbon (index 2 in sugar_match)
    phosphate_found = False
    phosphate_pattern = Chem.MolFromSmarts("[O][P](=O)([O])[O]")
    
    for sugar_match in sugar_matches:
        c3 = sugar_match[2]  # Third carbon (3' position)
        c5 = sugar_match[0]  # First carbon (5' position)
        
        # Check 5' position
        for bond in mol.GetAtomWithIdx(c5).GetBonds():
            if bond.GetBondType() == Chem.BondType.SINGLE:
                neighbor = bond.GetOtherAtomIdx(c5)
                if mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 15:  # Phosphorus
                    phosphate_found = True
                    break
        # Check 3' position
        for bond in mol.GetAtomWithIdx(c3).GetBonds():
            if bond.GetBondType() == Chem.BondType.SINGLE:
                neighbor = bond.GetOtherAtomIdx(c3)
                if mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 15:  # Phosphorus
                    phosphate_found = True
                    break
        if phosphate_found:
            break

    if not phosphate_found:
        return False, "No phosphate attached to 3' or 5' position"

    # 4. Exclude CoA-like structures with bridge phosphates
    # Look for thioester linkages (common in CoA derivatives)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("SCCNC(=O)")):
        return False, "Contains CoA-like thioester group"

    return True, "Contains nucleoside with phosphate at 3' or 5' position"