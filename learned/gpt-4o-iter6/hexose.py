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
    
    # Aldohexose pattern: check for aldehyde group (R-CHO)
    aldohexose_pattern = Chem.MolFromSmarts("[CX3H1](=O)[C]")
    if mol.HasSubstructMatch(aldohexose_pattern):
        return True, "Contains aldehyde group indicating aldohexose"
    
    # Ketohexose pattern: check for ketone group (R-CO-R)
    ketohexose_pattern = Chem.MolFromSmarts("[CX3](=O)[CX4]")
    if mol.HasSubstructMatch(ketohexose_pattern):
        return True, "Contains ketone group indicating ketohexose"

    # Cyclic furanose (5-membered) and pyranose (6-membered) patterns
    pyranose_pattern = Chem.MolFromSmarts("C1OC(CO)C(O)C(O)C1")
    furanose_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C1")
    
    if mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(furanose_pattern):
        return True, "Contains cyclic form indicating hexose"

    return False, "Does not match hexose structure criteria"