"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_long_chain_fatty_acid(smiles: str):
    """
    Determine if a molecule is a long-chain fatty acid based on its SMILES string.
    A long-chain fatty acid has a chain length ranging from C13 to C22 with a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Count total number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check if the number of carbons fits the long-chain fatty acid range
    if c_count < 13 or c_count > 22:
        return False, f"Carbon count out of range for long-chain fatty acid: {c_count} carbons"
    
    # Detect the longest carbon chain
    longest_chain_length = max(len(path) for path in rdmolops.FindAllPathsOfLengthN(mol, n=c_count, onlyBonds=False))
    
    # Validate the longest chain is within specified range
    if longest_chain_length < 13 or longest_chain_length > 22:
        return False, f"Longest carbon chain length out of range: {longest_chain_length} carbons"
    
    return True, "Contains carboxylic acid group and valid carbon chain length for long-chain fatty acid"