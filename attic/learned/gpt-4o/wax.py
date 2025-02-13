"""
Classifies: CHEBI:73702 wax
"""
from rdkit import Chem

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    A wax is typically an ester formed from long-chain fatty acids and long-chain alcohols.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester groups (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester group found"
    if len(ester_matches) > 2:
        return False, f"Too many ester groups, found {len(ester_matches)}"

    # Function to determine the longest carbon chain
    def get_longest_chain_length(atom, visited):
        visited.add(atom.GetIdx())
        chain_lengths = [
            get_longest_chain_length(neighbor, visited.copy()) + 1
            for neighbor in atom.GetNeighbors()
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited
        ]
        return max(chain_lengths, default=0)

    chain_lengths = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            chain_length = get_longest_chain_length(atom, set())
            chain_lengths.append(chain_length)

    # Ensure there are at least two long carbon chains
    long_chain_threshold = 12  # Needs adjustment based on specific requirements for waxes
    long_chains = [length for length in chain_lengths if length >= long_chain_threshold]

    if len(long_chains) < 2:
        return False, "Insufficient long carbon chains; need at least two"

    return True, "Contains long-chain molecules with one or two ester linkages characteristic of waxes"