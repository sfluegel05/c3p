"""
Classifies: CHEBI:23424 cyanides
"""
from rdkit import Chem

def is_cyanides(smiles: str):
    """
    Determines if a molecule is a cyanide (Salts and C-organyl derivatives of hydrogen cyanide, HC#N).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyanide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the cyano group (C#N)
    cyano_group = Chem.MolFromSmarts('C#N')
    if mol.HasSubstructMatch(cyano_group):
        return True, "Molecule contains a cyano group (C#N)"
    else:
        return False, "Molecule does not contain a cyano group (C#N)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23424',
                          'name': 'cyanides',
                          'definition': 'Salts and C-organyl derivatives of '
                                        'hydrogen cyanide, HC#N.',
                          'parents': ['CHEBI:35352']},
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
    'num_true_positives': 66,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}