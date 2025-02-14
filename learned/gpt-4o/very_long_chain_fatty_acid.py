"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem

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

    # Check for carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Identify the longest chain terminating at a carboxylic acid
    longest_chain_length = 0
    for atom in mol.GetAtoms():
        # Consider only carbon atom chains
        if atom.GetAtomicNum() == 6:
            chain_length = len(Chem.rdmolops.FindAtomEnvironmentOfRadiusN(mol, atom, 1000))
            if chain_length > longest_chain_length:
                longest_chain_length = chain_length

    # Check chain length criteria for very long-chain fatty acids (greater than 22 carbons)
    if longest_chain_length > 22:
        return True, f"Longest carbon chain length is {longest_chain_length}, which is greater than 22"

    return False, f"Longest carbon chain length is {longest_chain_length}, not greater than 22"