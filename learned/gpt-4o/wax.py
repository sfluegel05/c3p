"""
Classifies: CHEBI:73702 wax
"""
from rdkit import Chem

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    A wax is typically an ester formed between long-chain fatty acids and long-chain alcohols.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester group (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"

    # Use a search algorithm to determine the longest continuous carbon chain
    def get_longest_chain_length(atom, visited):
        # Depth-first search to find the longest path
        visited.add(atom.GetIdx())
        chain_lengths = [
            get_longest_chain_length(neighbor, visited.copy()) + 1
            for neighbor in atom.GetNeighbors()
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited
        ]
        return max(chain_lengths, default=0)

    max_chain_length = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            chain_length = get_longest_chain_length(atom, set())
            if chain_length > max_chain_length:
                max_chain_length = chain_length

    # Wax should have at least two long carbon chains
    long_chain_threshold = 12  # This could be adjusted as needed
    if max_chain_length < long_chain_threshold:
        return False, f"Insufficient long carbon chains, longest found is {max_chain_length}"

    # Check for exactly one ester linkage
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Incorrect number of ester groups, found {len(ester_matches)}"

    return True, "Contains long-chain molecules with one ester linkage characteristic of waxes"