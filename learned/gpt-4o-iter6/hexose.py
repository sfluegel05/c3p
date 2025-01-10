"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    Hexoses are six-carbon monosaccharides which may exist as aldohexoses or ketohexoses 
    in linear or cyclic forms.
  
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure there are exactly 6 carbon atoms in the hexose backbone
    backbone_pattern = Chem.MolFromSmarts("C-C-C-C-C-C")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "Does not contain a six-carbon backbone typical of hexoses"

    # Aldohexose pattern: check for aldehyde group at the first carbon of the backbone
    aldohexose_pattern = Chem.MolFromSmarts("C-C-C-C-C-C(=O)C")
    if mol.HasSubstructMatch(aldohexose_pattern):
        return True, "Contains aldehyde group indicating aldohexose"
    
    # Ketohexose pattern: check for ketone group at the second carbon of the backbone
    ketohexose_pattern = Chem.MolFromSmarts("C-C(=O)C-C-C-C")
    if mol.HasSubstructMatch(ketohexose_pattern):
        return True, "Contains ketone group indicating ketohexose"

    # Cyclic pyranose (6-membered) and furanose (5-membered) patterns need to consider more stereochemistry
    pyranose_pattern = Chem.MolFromSmarts("C1[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)[C@H]1")
    furanose_pattern = Chem.MolFromSmarts("C1[C@H](O)[C@H](O)[C@@H](O)[C@H]1")

    if mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(furanose_pattern):
        return True, "Contains cyclic form indicating hexose"

    return False, "Does not match hexose structure criteria"