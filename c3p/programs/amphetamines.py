"""
Classifies: CHEBI:35338 amphetamines
"""
from rdkit import Chem

def is_amphetamines(smiles: str):
    """
    Determines if a molecule is an amphetamine (based on the structure of the parent amphetamine 1-phenylpropan-2-amine).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amphetamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the 1-phenylpropan-2-amine core structure
    amphetamine_core = Chem.MolFromSmarts("CC(CC1=CC=CC=C1)N")
    if mol.HasSubstructMatch(amphetamine_core):
        return True, "Contains the 1-phenylpropan-2-amine core structure"
    
    return False, "Does not contain the 1-phenylpropan-2-amine core structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35338',
                          'name': 'amphetamines',
                          'definition': 'Amines that constitute a class of '
                                        'central nervous system stimulants '
                                        'based on the structure of the parent '
                                        'amphetamine 1-phenylpropan-2-amine.',
                          'parents': ['CHEBI:32952']},
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
    'num_true_negatives': 10,
    'num_false_negatives': 10,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}