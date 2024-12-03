"""
Classifies: CHEBI:22750 benzylisoquinoline alkaloid
"""
from rdkit import Chem

def is_benzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a benzylisoquinoline alkaloid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a benzylisoquinoline alkaloid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of isoquinoline ring system
    isoquinoline = Chem.MolFromSmarts('c1ccc2ncccc2c1')
    if not mol.HasSubstructMatch(isoquinoline):
        return False, "No isoquinoline ring system found"

    # Check for the presence of benzyl group
    benzyl = Chem.MolFromSmarts('c1ccccc1CC')
    if not mol.HasSubstructMatch(benzyl):
        return False, "No benzyl group found"

    # Check for the connection between benzyl group and isoquinoline system
    benzyl_isoquinoline = Chem.MolFromSmarts('c1ccc2ncccc2c1CC')
    if not mol.HasSubstructMatch(benzyl_isoquinoline):
        return False, "No benzylisoquinoline skeleton found"

    return True, "Molecule is a benzylisoquinoline alkaloid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22750',
                          'name': 'benzylisoquinoline alkaloid',
                          'definition': 'Any isoquinoline alkaloid based on a '
                                        'benzylisoquinoline skeleton.',
                          'parents': ['CHEBI:24921']},
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
    'num_true_negatives': 12,
    'num_false_negatives': 12,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}