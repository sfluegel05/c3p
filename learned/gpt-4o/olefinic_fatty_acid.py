"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
from rdkit import Chem

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid based on its SMILES string.
    An olefinic fatty acid contains at least one C=C double bond and a carboxylic acid group,
    or is part of an esterified or complex lipid molecule with olefinic characteristics.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an olefinic fatty acid or derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Search for C=C double bond (carbon backbone double bond typical of fatty acids)
    double_bond_pattern = Chem.MolFromSmarts("C=CC")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No C=C double bond in likely carbon backbone position found"
    
    # Broaden search for a long enough carbon chain
    # A pattern implying a flexible variable-length chain can help
    long_chain_pattern = Chem.MolFromSmarts("C(CC)C")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No suitable long carbon chain typical of fatty acids found"
    
    # Search for carboxylic acid group or esterification in variable contexts
    carboxylic_acid_or_ester_pattern = Chem.MolFromSmarts("C(=O)[O,OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_or_ester_pattern):
        return False, "No carboxylic acid or esterified group found in structure"

    return True, "Contains olefinic characteristics: C=C double bond, long carbon chain, and carboxylic acid or esterified group"