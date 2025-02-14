"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
"""
Classifies: CHEBI:50786 tertiary amine oxide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    A tertiary amine oxide is an N-oxide where there are three organic groups bonded to the nitrogen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine oxide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for N-oxide pattern ([N+](=O)[!O;!H])
    n_oxide_pattern = Chem.MolFromSmarts("[N+](=O)[!O;!H]")
    n_oxide_matches = mol.GetSubstructMatches(n_oxide_pattern)
    if len(n_oxide_matches) != 1:
        return False, f"Found {len(n_oxide_matches)} N-oxide groups, need exactly one"

    # Check if N-oxide is tertiary (3 organic groups attached)
    n_oxide_idx = n_oxide_matches[0][0]
    n_oxide_atom = mol.GetAtomWithIdx(n_oxide_idx)
    organic_neighbors = [nbr for nbr in n_oxide_atom.GetNeighbors() if nbr.GetAtomicNum() != 8 and nbr.GetAtomicNum() != 1]
    if len(organic_neighbors) != 3:
        return False, f"N-oxide has {len(organic_neighbors)} organic neighbors, need exactly 3"

    # Check if organic groups are alkyl/aryl
    alkyl_aryl_pattern = Chem.MolFromSmarts("[C;A]")
    for nbr in organic_neighbors:
        if not mol.GetAtomWithIdx(nbr.GetIdx()).HasSubstructMatch(alkyl_aryl_pattern):
            return False, "One of the organic groups is not alkyl or aryl"

    # Additional checks (e.g., connectivity, rings, etc.)
    # ...

    return True, "Molecule is a tertiary amine oxide"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:50786',
        'name': 'tertiary amine oxide',
        'definition': 'An N-oxide where there are three organic groups bonded to the nitrogen atom.',
        'parents': ['CHEBI:35428', 'CHEBI:50781']
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
    'num_true_positives': 442,
    'num_false_positives': 10,
    'num_true_negatives': 182458,
    'num_false_negatives': 3,
    'num_negatives': None,
    'precision': 0.9782608695652174,
    'recall': 0.9932038834951456,
    'f1': 0.9856380457521599,
    'accuracy': 0.9998479811417816
}