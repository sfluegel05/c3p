"""
Classifies: CHEBI:51149 xanthones
"""
from rdkit import Chem

def is_xanthones(smiles: str):
    """
    Determines if a molecule is a xanthone (based on a xanthone skeleton).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the xanthone core structure
    xanthone_smarts = 'O=C1C=CC2=CC=CC=C2OC3=C1C=CC=C3'
    xanthone_core = Chem.MolFromSmarts(xanthone_smarts)
    
    if xanthone_core is None:
        return False, "Invalid xanthone SMARTS pattern"

    if mol.HasSubstructMatch(xanthone_core):
        return True, "Molecule contains xanthone core structure"
    else:
        return False, "Molecule does not contain xanthone core structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51149',
                          'name': 'xanthones',
                          'definition': 'Any member of the class of  xanthenes '
                                        'based on a xanthone skeleton.',
                          'parents': ['CHEBI:38835', 'CHEBI:3992']},
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
    'num_false_negatives': 52,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}