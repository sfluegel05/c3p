"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains a steroid backbone
    steroid_backbone = Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CCC4C3(CCCC4)C')
    if not mol.HasSubstructMatch(steroid_backbone):
        return False, "No steroid backbone found"

    # Check for 3alpha-hydroxy group
    hydroxy_3alpha = Chem.MolFromSmarts('[C@H](O)[C@]1(CC[C@H]2[C@@H]3[C@H](O)C[C@@H]4C[C@H](O)CC[C@]4(C)[C@H]3C[C@H](O)[C@]12C)')
    if not mol.HasSubstructMatch(hydroxy_3alpha):
        return False, "No 3alpha-hydroxy group found"

    return True, "3alpha-hydroxy steroid found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36835',
                          'name': '3alpha-hydroxy steroid',
                          'definition': 'A 3-hydroxy steroid in which the '
                                        '3-hydroxy substituent is in the '
                                        'alpha-position.',
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
    'num_true_negatives': 13,
    'num_false_negatives': 13,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}