"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for 3beta-hydroxy steroid
    pattern = Chem.MolFromSmarts("[C@@H]1([C@H](O)C[C@@H]2CC[C@]3(C)C4CCC5=CC(=O)CCC5(C)[C@]4(C)CC[C@]3(C)[C@]2(C)[C@H]1C)")

    if mol.HasSubstructMatch(pattern):
        return True, "Molecule matches the 3beta-hydroxy steroid pattern"
    else:
        return False, "Molecule does not match the 3beta-hydroxy steroid pattern"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36836',
                          'name': '3beta-hydroxy steroid',
                          'definition': 'A 3-hydroxy steroid in which the '
                                        '3-hydroxy substituent is in the  '
                                        'beta-position.',
                          'parents': ['CHEBI:35681', 'CHEBI:36834']},
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
    'num_false_negatives': 38,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}