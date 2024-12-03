"""
Classifies: CHEBI:46955 anthracenes
"""
from rdkit import Chem

def is_anthracenes(smiles: str):
    """
    Determines if a molecule contains an anthracene skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains an anthracene skeleton, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the anthracene skeleton SMARTS pattern
    anthracene_smarts = "c1ccc2cc3ccccc3cc2c1"
    anthracene_pattern = Chem.MolFromSmarts(anthracene_smarts)

    if mol.HasSubstructMatch(anthracene_pattern):
        return True, "Contains anthracene skeleton"
    else:
        return False, "Does not contain anthracene skeleton"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:46955',
                          'name': 'anthracenes',
                          'definition': 'Compounds containing an anthracene '
                                        'skeleton.',
                          'parents': ['CHEBI:38032', 'CHEBI:51269']},
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
    'num_true_positives': 2,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 82,
    'precision': 0.6666666666666666,
    'recall': 0.023809523809523808,
    'f1': 0.04597701149425287,
    'accuracy': None}