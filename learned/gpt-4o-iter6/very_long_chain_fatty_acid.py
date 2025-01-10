"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determine if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid has a chain length greater than C22. 

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise
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
    
    # Find all carbon chains
    carbon_pattern = Chem.MolFromSmarts('[#6]')
    substruct_matches = mol.GetSubstructMatches(carbon_pattern)
    
    # Track the longest continuous chain of carbon atoms
    max_chain_length = 0
    
    for atom_tuple in substruct_matches:
        chain = [mol.GetAtomWithIdx(idx) for idx in atom_tuple]
        
        if all(atom.GetSymbol() == 'C' for atom in chain):
            num_carbons_in_chain = len(chain)
            max_chain_length = max(max_chain_length, num_carbons_in_chain)
    
    # Determine if the longest carbon chain is a very long-chain fatty acid
    if max_chain_length > 22:
        return True, f"Contains {max_chain_length} carbons in a continuous chain, qualifies as a very long-chain fatty acid"
    else:
        return False, f"Contains {max_chain_length} carbons in longest chain, does not exceed C22"

# Example usage:
# result, reason = is_very_long_chain_fatty_acid("C(O)(=O)CCCCCCCCCCCCCCCCCCCCCC(O)=O")
# print(result, reason)