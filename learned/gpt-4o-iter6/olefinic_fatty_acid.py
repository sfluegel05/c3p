"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
from rdkit import Chem

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid.
    An olefinic fatty acid is characterized by having a carboxylic acid group
    and at least one carbon-carbon double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an olefinic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Detect carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Detect carbon-carbon double bond (C=C)
    cc_double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(cc_double_bond_pattern):
        return False, "No carbon-carbon double bond found"

    return True, "Contains both a carboxylic acid group and at least one C=C double bond"