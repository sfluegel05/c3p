"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem

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
    
    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Use RDKit's built-in method to find the longest contiguous carbon chain
    # Use getLongestAliphaticChain interpret long chains in fatty acids well
    path = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    atom_weights = []
    while path:
        next_longest_cpath = Chem.GetLongestAliphaticChain(mol)
        atom_weights.append(next_longest_cpath)
    max_chain_length = max(atom_weights)

    # Determine if the longest carbon chain is a very long-chain fatty acid
    if max_chain_length > 22:
        return True, f"Contains {max_chain_length} carbons in a continuous chain, qualifies as a very long-chain fatty acid"
    else:
        return False, f"Contains {max_chain_length} carbons in longest chain, does not exceed C22"

# Example usage:
# result, reason = is_very_long_chain_fatty_acid("C(O)(=O)CCCCCCCCCCCCCCCCCCCCCC(O)=O")
# print(result, reason)