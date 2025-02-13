"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: CHEBI:27388 very long-chain fatty acid
"""
from rdkit import Chem

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid is defined as a fatty acid with a chain length greater than C22.
    Ultra-long-chain fatty acids have a chain length greater than C27.

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

    # Look for carboxylic acid group (-COOH)
    carboxylic_acid = Chem.MolFromSmarts('C(=O)[OH]')
    matches = mol.GetSubstructMatches(carboxylic_acid)
    if len(matches) == 0:
        return False, "No carboxylic acid group found"

    # Identify all possible chains starting from carboxyl carbon
    max_chain_length = 0
    for match in matches:
        carboxyl_carbon_idx = match[0]
        visited = set()

        def dfs(atom_idx, visited):
            visited = visited.copy()
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)

            length = 1  # Current carbon
            for neighbor in atom.GetNeighbors():
                nbr_idx = neighbor.GetIdx()
                neighbor_atom = mol.GetAtomWithIdx(nbr_idx)
                # Traverse to neighboring carbons not yet visited
                if neighbor_atom.GetAtomicNum() == 6 and nbr_idx not in visited:
                    length = max(length, 1 + dfs(nbr_idx, visited))
            return length

        chain_length = dfs(carboxyl_carbon_idx, visited)
        if chain_length > max_chain_length:
            max_chain_length = chain_length

    if max_chain_length == 0:
        return False, "No suitable carbon chain found"

    # Check if chain length exceeds 22 carbons
    if max_chain_length > 22:
        if max_chain_length > 27:
            return True, f"Ultra-long-chain fatty acid with chain length C{max_chain_length}"
        else:
            return True, f"Very long-chain fatty acid with chain length C{max_chain_length}"
    else:
        return False, f"Chain length is C{max_chain_length}, not greater than C22"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:27388',
                             'name': 'very long-chain fatty acid',
                             'definition': 'A fatty acid which has a chain length greater than C22. Very long-chain fatty acids which have a chain length greater than C27 are also known as ultra-long-chain fatty acids.',
                             'parents': []},
       'config': {   'llm_model_name': 'lbl/claude-sonnet',
                     'f1_threshold': 0.8,
                     'max_attempts': 5,
                     'max_positive_instances': None,
                     'max_positive_to_test': None,
                     'max_negative_to_test': None,
                     'max_positive_in_prompt': 50,
                     'max_negative_in_prompt': 20,
                     'max_instances_in_prompt': 100,
                     'test_proportion': 0.1},
       'message': None,
       'attempt': 3,
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
       'accuracy': None}