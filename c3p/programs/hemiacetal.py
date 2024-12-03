"""
Classifies: CHEBI:5653 hemiacetal
"""
from rdkit import Chem

def is_hemiacetal(smiles: str):
    """
    Determines if a molecule is a hemiacetal.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hemiacetal, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a hemiacetal
    hemiacetal_pattern = Chem.MolFromSmarts('[#6][#6](O)[#6][O]')

    if mol.HasSubstructMatch(hemiacetal_pattern):
        return True, "Molecule contains a hemiacetal group"
    else:
        return False, "Molecule does not contain a hemiacetal group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:5653',
                          'name': 'hemiacetal',
                          'definition': 'A compound having the general formula '
                                        "RR'C(OH)OR'' (R'' =/= H).",
                          'parents': ['CHEBI:30879']},
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
    'num_true_positives': 9,
    'num_false_positives': 2,
    'num_true_negatives': 12,
    'num_false_negatives': 5,
    'precision': 0.8181818181818182,
    'recall': 0.6428571428571429,
    'f1': 0.7200000000000001,
    'accuracy': None}