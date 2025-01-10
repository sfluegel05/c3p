"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
from rdkit import Chem

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid based on its SMILES string.
    A long-chain fatty acid is a fatty acid with a chain length ranging from C13 to C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxylic acid groups (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-1]")
    carboxyl_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"

    # Function to find the longest carbon chain from a given atom
    def get_longest_chain_length(atom_idx, visited):
        """
        Recursively finds the longest carbon chain length from the given atom.

        Args:
            atom_idx (int): Index of the current atom
            visited (set): Set of visited atom indices

        Returns:
            int: Longest chain length from the current atom
        """
        visited.add(atom_idx)
        max_length = 1  # Count current atom

        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in visited and neighbor.GetAtomicNum() == 6:
                chain_length = 1 + get_longest_chain_length(neighbor_idx, visited)
                if chain_length > max_length:
                    max_length = chain_length
        visited.remove(atom_idx)
        return max_length

    # Find the longest carbon chain attached to any carboxyl carbon
    max_chain_length = 0
    for match in carboxyl_matches:
        carboxyl_carbon_idx = match[0]  # First atom in the pattern is the carboxyl carbon
        visited = set()
        chain_length = get_longest_chain_length(carboxyl_carbon_idx, visited)
        if chain_length > max_chain_length:
            max_chain_length = chain_length

    # Check if the chain length is within the specified range
    if 13 <= max_chain_length <= 22:
        return True, f"Longest carbon chain length is {max_chain_length}"
    else:
        return False, f"Longest carbon chain length is {max_chain_length}, which is outside the range 13-22"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:15904',
        'name': 'long-chain fatty acid',
        'definition': 'A fatty acid with a chain length ranging from C13 to C22.',
        'parents': ['CHEBI:15718']  # fatty acid
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}