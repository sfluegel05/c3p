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
    
    # Search for C=C double bond in a long enough alkyl chain typical of fatty acids
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 1:
        return False, "No C=C double bond found"
    
    # Verify the presence of a carbon chain
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCC")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long carbon chain found"
    
    # Search for carboxylic acid group or esterified form
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    esterified_pattern = Chem.MolFromSmarts("COC(=O)C")
    
    if not (mol.HasSubstructMatch(carboxylic_acid_pattern) or 
            mol.HasSubstructMatch(esterified_pattern)):
        return False, "No carboxylic acid or esterified group found"
    
    return True, "Contains at least one C=C double bond, long carbon chain, and carboxylic acid or esterified group"