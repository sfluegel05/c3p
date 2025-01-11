"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
"""
Classifies: CHEBI:59507 olefinic fatty acid
"""
from rdkit import Chem

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid based on its SMILES string.
    An olefinic fatty acid is a fatty acid containing at least one C=C double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an olefinic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)OH or -C(=O)O-)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for at least one C=C double bond
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No C=C double bond found"

    return True, "Contains a carboxylic acid group and at least one C=C double bond"