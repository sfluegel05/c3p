"""
Classifies: CHEBI:23238 chromones
"""
from rdkit import Chem

def is_chromones(smiles: str):
    """
    Determines if a molecule is a chromone (1,4-benzopyrone skeleton and its substituted derivatives).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chromone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the 1,4-benzopyrone skeleton
    chromone_pattern = Chem.MolFromSmarts('O=C1C=CC2=CC=CC=C2O1')
    if chromone_pattern is None:
        return False, "Invalid SMARTS pattern for chromone"

    # Check if the molecule matches the chromone pattern
    if mol.HasSubstructMatch(chromone_pattern):
        return True, "Molecule contains the 1,4-benzopyrone skeleton"
    else:
        return False, "Molecule does not contain the 1,4-benzopyrone skeleton"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23238',
                          'name': 'chromones',
                          'definition': 'A chromenone that consists of a '
                                        '1,4-benzopyrone skeleton and its '
                                        'substituted derivatives thereof.',
                          'parents': ['CHEBI:38445']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 13-14: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}