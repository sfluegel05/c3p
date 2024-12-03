"""
Classifies: CHEBI:36786 tetralins
"""
from rdkit import Chem

def is_tetralins(smiles: str):
    """
    Determines if a molecule contains a tetralin skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a tetralin skeleton, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the tetralin substructure as a SMARTS pattern
    tetralin_smarts = "C1CCC2=CC=CC=C2C1"
    tetralin_substructure = Chem.MolFromSmarts(tetralin_smarts)

    if mol.HasSubstructMatch(tetralin_substructure):
        return True, "Molecule contains a tetralin skeleton"
    else:
        return False, "Molecule does not contain a tetralin skeleton"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36786',
                          'name': 'tetralins',
                          'definition': 'Compounds containing a tetralin '
                                        'skeleton.',
                          'parents': ['CHEBI:33836', 'CHEBI:36785']},
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 36,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}