"""
Classifies: CHEBI:38672 flavans
"""
from rdkit import Chem

def is_flavans(smiles: str):
    """
    Determines if a molecule is a flavan (3,4-dihydro-2-aryl-2H-1-benzopyran skeleton and its substituted derivatives).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavan, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavan core structure
    flavan_core = Chem.MolFromSmarts('C1C2=CC=CC=C2C(C(O)C3=CC=CC=C3)O1')

    if flavan_core is None:
        return False, "Invalid flavan core structure"

    # Check if molecule contains the flavan core
    if mol.HasSubstructMatch(flavan_core):
        return True, "Molecule contains the flavan core structure"
    else:
        return False, "Molecule does not contain the flavan core structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38672',
                          'name': 'flavans',
                          'definition': 'Any flavonoid with a '
                                        '3,4-dihydro-2-aryl-2H-1-benzopyran '
                                        'skeleton and its substituted '
                                        'derivatives.',
                          'parents': [   'CHEBI:26004',
                                         'CHEBI:38443',
                                         'CHEBI:47916']},
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
    'num_false_negatives': 68,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}