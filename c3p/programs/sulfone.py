"""
Classifies: CHEBI:35850 sulfone
"""
from rdkit import Chem

def is_sulfone(smiles: str):
    """
    Determines if a molecule is a sulfone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define sulfone substructure
    sulfone_smarts = '[#6][S](=O)(=O)[#6]'
    sulfone = Chem.MolFromSmarts(sulfone_smarts)

    if mol.HasSubstructMatch(sulfone):
        return True, "Molecule contains sulfone group"
    else:
        return False, "Molecule does not contain sulfone group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35850',
                          'name': 'sulfone',
                          'definition': 'An organosulfur compound having the '
                                        'structure RS(=O)2R (R =/= H).',
                          'parents': ['CHEBI:33261']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 30-31: malformed \\N character escape (<string>, line '
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