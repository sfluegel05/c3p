"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determine if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid has a chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group presence
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Finding and analyzing chains: longest single, after excluding specific groups not part of the chain.
    chains = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetHybridization() < Chem.rdchem.HybridizationType.SP:  # Only sp3 or less
            bfs_path = rdmolops.GetShortestPath(mol, 0, atom.GetIdx())
            chain_len = sum(1 for _ in bfs_path) - 1  # excluding the acid C(=O)
            chains.append(chain_len)
    
    # Get maximum chain length
    if chains:
        longest_chain = max(chains)
    else:
        return False, "No carbon chains found"

    # Determine if the longest carbon chain qualifies
    if longest_chain > 22:
        return True, f"Contains {longest_chain} carbons in a continuous chain, qualifies as a very long-chain fatty acid"
    else:
        return False, f"Contains {longest_chain} carbons in the longest chain, does not exceed C22"