"""
Classifies: CHEBI:46700 benzoxazole
"""
from rdkit import Chem

def is_benzoxazole(smiles: str):
    """
    Determines if a molecule is a benzoxazole (compounds based on a fused 1,2- or 1,3-oxazole and benzene bicyclic ring skeleton).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a benzoxazole, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for benzoxazole
    benzoxazole_smarts = "[cH]1[cH][cH][cH][cH][cH]1[c]2[o,n][c][n,o]2"
    
    # Check if the molecule matches the benzoxazole pattern
    benzoxazole_pattern = Chem.MolFromSmarts(benzoxazole_smarts)
    if mol.HasSubstructMatch(benzoxazole_pattern):
        return True, "Molecule is a benzoxazole"
    else:
        return False, "Molecule does not match benzoxazole pattern"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:46700',
                          'name': 'benzoxazole',
                          'definition': 'Compounds based on a fused 1,2- or '
                                        '1,3-oxazole and benzene bicyclic ring '
                                        'skeleton.',
                          'parents': [   'CHEBI:27171',
                                         'CHEBI:38101',
                                         'CHEBI:38104']},
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