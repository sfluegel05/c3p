"""
Classifies: CHEBI:39447 pyrimidines
"""
from rdkit import Chem

def is_pyrimidines(smiles: str):
    """
    Determines if a molecule has a pyrimidine as part of its structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a pyrimidine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the pyrimidine substructure
    pyrimidine_smarts = "c1cncnc1"
    pyrimidine = Chem.MolFromSmarts(pyrimidine_smarts)

    if mol.HasSubstructMatch(pyrimidine):
        return True, "Pyrimidine structure found"
    else:
        return False, "No pyrimidine structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:39447',
                          'name': 'pyrimidines',
                          'definition': 'Any compound having a pyrimidine as '
                                        'part of its structure.',
                          'parents': ['CHEBI:38313']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 14-15: malformed \\N character escape (<string>, line '
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