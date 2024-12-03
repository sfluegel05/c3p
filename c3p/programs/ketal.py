"""
Classifies: CHEBI:59777 ketal
"""
from rdkit import Chem

def is_ketal(smiles: str):
    """
    Determines if a molecule is a ketal.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ketal, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a ketal
    ketal_pattern = Chem.MolFromSmarts("[C](O)(O)[C]")

    if mol.HasSubstructMatch(ketal_pattern):
        return True, "Molecule contains a ketal group"
    else:
        return False, "Molecule does not contain a ketal group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:59777',
                          'name': 'ketal',
                          'definition': 'An acetal of formula R2C(OR)2 (R =/= '
                                        'H) derived from a ketone by '
                                        'replacement of the oxo group by two '
                                        'hydrocarbyloxy groups. The class name '
                                        "'ketals', once abandoned by IUPAC, "
                                        'has been reinstated as a subclass of '
                                        'acetals.',
                          'parents': ['CHEBI:59769']},
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
    'num_true_positives': 33,
    'num_false_positives': 7,
    'num_true_negatives': 13,
    'num_false_negatives': 0,
    'precision': 0.825,
    'recall': 1.0,
    'f1': 0.9041095890410958,
    'accuracy': None}