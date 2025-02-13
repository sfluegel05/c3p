"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid has a carbon chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure the molecule has a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    carboxyl_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"
    
    max_chain_length = 0
    
    # Measure chain length from each carboxylic acid group found
    for index, match in enumerate(carboxyl_matches):
        # Get the carbon atom in the carboxylic acid (usually the first in match)
        carbon_in_carboxyl = match[0]
        
        # Identify the longest chain starting from carbon in the carboxylic group
        chain_length = rdmolops.GetShortestPathLength(mol, carbon_in_carboxyl)
        
        if chain_length > max_chain_length:
            max_chain_length = chain_length
    
    # Since carboxylic acid carbon itself is included in the length, reduce by 1
    max_chain_length -= 1

    # Check chain length criteria for very long-chain fatty acids (greater than 22 carbons)
    if max_chain_length > 22:
        return True, f"Longest carbon chain length is {max_chain_length}, which is greater than 22"

    return False, f"Longest carbon chain length is {max_chain_length}, not greater than 22"