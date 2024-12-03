"""
Classifies: CHEBI:17478 aldehyde
"""
from rdkit import Chem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldehyde, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the aldehyde group
    aldehyde_pattern = Chem.MolFromSmarts('[CX3H1](=O)[#6]')
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Aldehyde group detected"
    else:
        return False, "Aldehyde group not detected"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17478',
                          'name': 'aldehyde',
                          'definition': 'A compound RC(=O)H, in which a '
                                        'carbonyl group is bonded to one '
                                        'hydrogen atom and to one R group.',
                          'parents': ['CHEBI:36586']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[13:21:59] WARNING: not removing hydrogen atom without '
             'neighbors\n',
    'stdout': '',
    'num_true_positives': 102,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}