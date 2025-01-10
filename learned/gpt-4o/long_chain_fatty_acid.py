"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
from rdkit import Chem

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
    
    # Identify the longest carbon chain
    chains = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.SINGLE:
            path = Chem.rdmolops.GetShortestPath(mol, bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            chain_length = sum(1 for atom_idx in path if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6)
            chains.append(chain_length)
    
    # Get the maximum chain length found
    longest_chain_length = max(chains) if chains else 0
    
    # Check if the longest chain length fits the long-chain fatty acid range
    if longest_chain_length < 13 or longest_chain_length > 22:
        return False, f"Longest carbon chain length out of range: {longest_chain_length} carbons"
    
    return True, "Contains carboxylic acid group and valid carbon chain length for long-chain fatty acid"