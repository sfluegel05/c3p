"""
Classifies: CHEBI:37912 hydroxycoumarin
"""
from rdkit import Chem

def is_hydroxycoumarin(smiles: str):
    """
    Determines if a molecule is a hydroxycoumarin (Any coumarin carrying at least one hydroxy substituent).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxycoumarin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the coumarin core (benzopyran-2-one)
    coumarin_smarts = '[O]=C1OC2=CC=CC=C2C=C1'
    coumarin = Chem.MolFromSmarts(coumarin_smarts)
    
    if not mol.HasSubstructMatch(coumarin):
        return False, "No coumarin core structure found"

    # Check for hydroxy substituent
    hydroxy_smarts = '[OH]'
    hydroxy = Chem.MolFromSmarts(hydroxy_smarts)
    
    if mol.HasSubstructMatch(hydroxy):
        return True, "Hydroxycoumarin with at least one hydroxy substituent"
    else:
        return False, "No hydroxy substituent found on coumarin core"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37912',
                          'name': 'hydroxycoumarin',
                          'definition': 'Any coumarin carrying at least one '
                                        'hydroxy substituent.',
                          'parents': ['CHEBI:23403']},
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
    'num_true_negatives': 19,
    'num_false_negatives': 19,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}