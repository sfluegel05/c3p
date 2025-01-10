"""
Classifies: CHEBI:35366 fatty acid
"""
"""
Classifies: CHEBI:28868 fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    A fatty acid is an aliphatic monocarboxylic acid with a chain of 4 to 28 carbons,
    which may be saturated or unsaturated, and can include hydroxyl groups or minor substitutions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxylic acid groups: [C](=O)[O-] or [C](=O)[OH]
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H1,OX1H0-]")
    carboxylic_acids = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acids) != 1:
        return False, f"Found {len(carboxylic_acids)} carboxylic acid groups, need exactly 1"

    # Get the carboxyl carbon atom index
    carboxyl_c_index = carboxylic_acids[0][0]
    visited = set()

    # Define atoms allowed in fatty acids
    allowed_atoms = {6, 1, 8, 16, 17, 35, 53}  # C, H, O, S, Cl, Br, I

    # Function to recursively find the longest carbon chain from a given atom
    def find_longest_chain(atom_idx, length):
        visited.add(atom_idx)
        max_length = length
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            n_atomic_num = neighbor.GetAtomicNum()
            bond = mol.GetBondBetweenAtoms(atom_idx, n_idx)
            # Allow traversal through C, O, S atoms (common in fatty acids)
            if n_atomic_num in allowed_atoms and n_idx not in visited:
                chain_length = find_longest_chain(n_idx, length + (1 if n_atomic_num == 6 else 0))
                if chain_length > max_length:
                    max_length = chain_length
        visited.remove(atom_idx)
        return max_length

    # Find the longest carbon chain starting from the carboxyl carbon
    chain_length = find_longest_chain(carboxyl_c_index, 0)

    # Check if chain length is within typical fatty acid range
    if chain_length < 4:
        return False, f"Chain length {chain_length} is too short for a fatty acid"
    elif chain_length > 28:
        return False, f"Chain length {chain_length} is too long for typical fatty acids"

    # Check if molecule contains only allowed elements (C, H, O, S, halogens)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, f"Contains disallowed element: {atom.GetSymbol()}"

    return True, f"Contains a single carboxylic acid group and an aliphatic chain of {chain_length} carbons"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:28868',
        'name': 'fatty acid',
        'definition': 'Any aliphatic monocarboxylic acid derived from or contained in esterified form in an animal or vegetable fat, oil or wax. Natural fatty acids commonly have a chain of 4 to 28 carbons (usually unbranched and even-numbered), which may be saturated or unsaturated. By extension, the term is sometimes used to embrace all acyclic aliphatic carboxylic acids.',
        'parents': ['CHEBI:25384', 'CHEBI:17436']
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
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199
}