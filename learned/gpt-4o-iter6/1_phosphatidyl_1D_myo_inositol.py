"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    The molecule should have a 1D-myo-inositol with a phosphatidyl group at position 1.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-phosphatidyl-1D-myo-inositol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for 1D-myo-inositol ring with appropriate stereochemistry
    inositol_smarts = "[C@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"
    inositol_pattern = Chem.MolFromSmarts(inositol_smarts)
    if inositol_pattern is None or not mol.HasSubstructMatch(inositol_pattern):
        return False, "No 1D-myo-inositol moiety found"

    # SMARTS for phosphate ester linkage; pattern revised for more flexibility
    phosphate_smarts = "C1OC([O,P]=O)=O"
    phosphate_linkage_pattern = Chem.MolFromSmarts(phosphate_smarts)
    if phosphate_linkage_pattern is None or not mol.HasSubstructMatch(phosphate_linkage_pattern):
        return False, "No phosphate ester linkage at position 1 of inositol"

    # Ensure the presence of glycerol backbone typically seen in phosphatidyl groups
    glycerol_smarts = "OCC(OP)=O"
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    if glycerol_pattern is None or not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone for phosphatidyl linkage"

    # Ensure phosphate group linkage to 1-position of inositol ring
    inositol_with_phosphate_smarts = "[C@H]1(O)[C@H](O)CC[P](=O)(O)O[C@H]1O"
    inositol_with_phosphate_pattern = Chem.MolFromSmarts(inositol_with_phosphate_smarts)
    if inositol_with_phosphate_pattern is None or not mol.HasSubstructMatch(inositol_with_phosphate_pattern):
        return False, "Phosphate group not correctly linked to inositol"

    return True, "Contains 1-phosphatidyl-1D-myo-inositol moiety"