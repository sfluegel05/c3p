"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    The molecule should have a 1D-myo-inositol with phosphatidyl group at position 1.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-phosphatidyl-1D-myo-inositol, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define 1D-myo-inositol SMARTS pattern: 6-membered ring with OH groups in specific stereo arrangement
    inositol_pattern = Chem.MolFromSmarts("[C@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No 1D-myo-inositol moiety found"

    # Define phosphatidylgroup SMARTS pattern (phosphate linkage attached to inositol)
    phosphatidyl_pattern = Chem.MolFromSmarts("O[C@@H]1(O[P](O)(=O)O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(phosphatidyl_pattern):
        return False, "No phosphatidyl group at position 1 of inositol"

    # Check glycerol backbone typically associated with phosphatidyls
    glycerol_pattern = Chem.MolFromSmarts("[C@H](COC(=O)[C,C]")  # A simple pattern indicating glycerol backbone
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    return True, "Contains 1-phosphatidyl-1D-myo-inositol moiety"