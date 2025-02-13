"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
"""
Classifies: CHEBI:51165 ultra-long-chain fatty acid
An ultra-long-chain fatty acid is any very long-chain fatty acid which has a chain length greater than C27.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ultra-long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Does not contain a carboxylic acid group"
    
    # Count carbon chain length
    carbon_chain_pattern = Chem.MolFromSmarts("[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]")
    chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    max_chain_length = 0
    for match in chain_matches:
        chain_length = len(match) + 1  # Add 1 for the carboxylic acid carbon
        max_chain_length = max(max_chain_length, chain_length)
    
    if max_chain_length <= 27:
        return False, f"Carbon chain length is only {max_chain_length}, need greater than 27"
    
    return True, f"Contains a carbon chain of length {max_chain_length} with a carboxylic acid group"