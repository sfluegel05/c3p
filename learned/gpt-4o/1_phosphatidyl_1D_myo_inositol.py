"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-phosphatidyl-1D-myo-inositol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")
    
    # Define the pattern for 1D-myo-inositol with correct stereochemistry
    inositol_pattern = Chem.MolFromSmarts("[C@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1(O)")
    if not mol.HasSubstructMatch(inositol_pattern):
        return (False, "1D-myo-inositol pattern not found or incorrect stereochemistry")

    # Define SMARTS pattern for phosphatidyl group with correct connectivity
    # Pattern is simplified: focus on glycerol backbone esterified with acyl chains and phosphate group
    phosphatidyl_pattern = Chem.MolFromSmarts("[C@H](CO[P](=O)(O)O)(O)COC(=O)C")
    if not mol.HasSubstructMatch(phosphatidyl_pattern):
        return (False, "Phosphatidyl group pattern not detected or incorrect")

    # Check the connection between phosphatidyl group and inositol
    inositol_match_atoms = mol.GetSubstructMatch(inositol_pattern)
    phosphatidyl_match_atoms = mol.GetSubstructMatch(phosphatidyl_pattern)

    if not inositol_match_atoms or not phosphatidyl_match_atoms:
        return (False, "Phosphatidyl or inositol pattern not fully matched in the molecule")

    # Verify the phosphatidyl connection at position 1 of inositol. 
    # Check if the connection involves position 1 carbon of the inositol mapped pattern
    inositol_atom_1 = inositol_match_atoms[0]  # Assuming SMARTS maps correctly
    phosphorus_atom = phosphatidyl_match_atoms[2]  # Phosphoryl group index in the phosphatidyl pattern
    
    for neighbor in mol.GetAtomWithIdx(inositol_atom_1).GetNeighbors():
        if neighbor.GetIdx() == phosphorus_atom:
            return (True, "Valid 1-phosphatidyl-1D-myo-inositol detected")
    
    return (False, "Phosphatidyl group not attached to position 1 of inositol")