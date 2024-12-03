"""
Classifies: CHEBI:24436 guanidines
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_guanidines(smiles: str):
    """
    Determines if a molecule is a guanidine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a guanidine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    guanidine_pattern = Chem.MolFromSmarts('NC(=N)N')
    if not mol.HasSubstructMatch(guanidine_pattern):
        return False, "No guanidine group found"

    return True, "Contains guanidine group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24436',
                          'name': 'guanidines',
                          'definition': 'Any organonitrogen compound '
                                        'containing a carbamimidamido '
                                        '(guanidino) group. Guanidines have '
                                        'the general structure '
                                        '(R(1)R(2)N)(R(3)R(4)N)C=N-R(5) and '
                                        'are related structurally to amidines '
                                        'and ureas.',
                          'parents': ['CHEBI:35352']},
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
    'num_true_positives': 18,
    'num_false_positives': 0,
    'num_true_negatives': 19,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9473684210526315,
    'f1': 0.972972972972973,
    'accuracy': None}