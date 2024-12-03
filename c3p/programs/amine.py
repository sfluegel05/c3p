"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for amines
    amine_patterns = [
        '[NX3;H2]',  # Primary amine
        '[NX3;H1][CX4]',  # Secondary amine
        '[NX3][CX4][CX4]',  # Tertiary amine
    ]

    for pattern in amine_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return True, "Molecule contains an amine group"

    return False, "Molecule does not contain an amine group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:32952',
                          'name': 'amine',
                          'definition': 'A compound formally derived from '
                                        'ammonia by replacing one, two or '
                                        'three hydrogen atoms by hydrocarbyl '
                                        'groups.',
                          'parents': ['CHEBI:50047']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 81,
    'num_false_positives': 13,
    'num_true_negatives': 7,
    'num_false_negatives': 5,
    'precision': 0.8617021276595744,
    'recall': 0.9418604651162791,
    'f1': 0.8999999999999999,
    'accuracy': None}