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

    # Ensure there are exactly 6 carbons
    carbon_pattern = Chem.MolFromSmarts("[#6]")
    c_count = len(mol.GetSubstructMatches(carbon_pattern))
    if c_count != 6:
        return False, f"Incorrect number of carbons for hexose: found {c_count}, expected exactly 6"

    # Aldohexose pattern: check for aldehyde group at position 1
    aldohexose_pattern = Chem.MolFromSmarts("[#6X3H1](=O)[#6]")
    if mol.HasSubstructMatch(aldohexose_pattern):
        return True, "Contains aldehyde group indicating aldohexose"
    
    # Ketohexose pattern: check for ketone group at position 2
    ketohexose_pattern = Chem.MolFromSmarts("[#6][#6X3](=O)[#6]")
    if mol.HasSubstructMatch(ketohexose_pattern):
        return True, "Contains ketone group indicating ketohexose"

    # Cyclic pyranose (6-membered) pattern with stereochemistry
    pyranose_pattern = Chem.MolFromSmarts("C1[C@H](O)[C@@H](O)C(O)[C@H](O)O1")
    # Cyclic furanose (5-membered) pattern with stereochemistry
    furanose_pattern = Chem.MolFromSmarts("C1[C@H](O)[C@@H](O)C(O)O1")
    
    if mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(furanose_pattern):
        return True, "Contains cyclic form indicating hexose"

    return False, "Does not match hexose structure criteria"