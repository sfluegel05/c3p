"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: epoxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid is a heterocyclic fatty acid containing an epoxide ring as part of its structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, 'Invalid SMILES string'

    # Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H1,H0-]')
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_matches:
        return False, "No carboxylic acid group found"

    # Assume the first carboxyl carbon as the starting point
    carboxyl_carbon_idx = carboxylic_matches[0][0]

    # Check for epoxide ring (three-membered cyclic ether)
    epoxide_pattern = Chem.MolFromSmarts('[C;R]1-[O;R]-[C;R]1')
    epoxide_matches = mol.GetSubstructMatches(epoxide_pattern)
    if not epoxide_matches:
        return False, "No epoxide ring found"

    # Collect epoxide carbon atom indices
    epoxide_carbon_indices = set()
    for match in epoxide_matches:
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:  # Carbon
                epoxide_carbon_indices.add(idx)

    # Function to find path from carboxyl carbon to epoxide carbons
    def find_path_to_epoxide(mol, start_idx, target_indices):
        visited = set()
        queue = [(start_idx, 0)]
        while queue:
            current_idx, path_length = queue.pop(0)
            if current_idx in target_indices:
                return path_length + 1  # Include current atom
            visited.add(current_idx)
            current_atom = mol.GetAtomWithIdx(current_idx)
            for neighbor in current_atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                neighbor_atom = mol.GetAtomWithIdx(neighbor_idx)
                if neighbor_idx not in visited and neighbor_atom.GetAtomicNum() == 6:
                    queue.append((neighbor_idx, path_length + 1))
        return None

    # Find path length
    path_length = find_path_to_epoxide(mol, carboxyl_carbon_idx, epoxide_carbon_indices)
    if path_length is None:
        return False, "No continuous carbon chain connecting carboxyl group and epoxide ring"

    # Check if the chain length is at least 8 carbons
    if path_length < 8:
        return False, f"Carbon chain length is {path_length}, which is less than 8 carbons"

    return True, "Molecule contains a carboxylic acid group and an epoxide ring within a hydrocarbon chain"


__metadata__ = {   'chemical_class': {   'id': None,
                                  'name': 'epoxy fatty acid',
                                  'definition': 'A heterocyclic fatty acid containing an epoxide ring as part of its structure.',
                                  'parents': ['fatty acid', 'heterocyclic compound']},
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
            'attempt': 1,
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