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
    
    # Look for a carboxylic acid group (C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Use Breadth First Search or Depth First Search method for longest carbon chain
    longest_chain_length = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:  # Skip non-carbon atoms
            continue
        chains = rdmolops.GetShortestPath(mol, atom.GetIdx())
        for chain in chains:
            chain_length = sum(1 for atom_idx in chain if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6)
            longest_chain_length = max(longest_chain_length, chain_length)
    
    # Check if the longest chain length fits the long-chain fatty acid range
    if longest_chain_length < 13 or longest_chain_length > 22:
        return False, f"Longest carbon chain length out of range: {longest_chain_length} carbons"
    
    return True, "Contains carboxylic acid group and valid carbon chain length for long-chain fatty acid"