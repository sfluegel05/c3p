"""
Classifies: CHEBI:35915 sterol ester
"""
from rdkit import Chem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester.
    
    A sterol ester is defined as a steroid ester obtained by formal condensation of the carboxy group of any carboxylic acid with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for ester functional group
    ester_found = False
    ester_smarts = Chem.MolFromSmarts("C(=O)O")
    if mol.HasSubstructMatch(ester_smarts):
        ester_found = True
    
    if not ester_found:
        return False, "No ester functional group found"
    
    # Check for sterol backbone
    sterol_found = False
    sterol_smarts = Chem.MolFromSmarts("[C@@H]1CC[C@H]2[C@@H]3CC=C4C[C@H](CC[C@]4(C)[C@H]3CC[C@]12C)")
    if mol.HasSubstructMatch(sterol_smarts):
        sterol_found = True
    
    if not sterol_found:
        return False, "No sterol backbone found"
    
    return True, "Sterol ester found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35915',
                          'name': 'sterol ester',
                          'definition': 'A steroid ester obtained by formal '
                                        'condensation of the carboxy group of '
                                        'any carboxylic acid with the '
                                        '3-hydroxy group of a sterol.',
                          'parents': ['CHEBI:33308', 'CHEBI:47880']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 7,
    'num_false_positives': 0,
    'num_true_negatives': 12,
    'num_false_negatives': 5,
    'precision': 1.0,
    'recall': 0.5833333333333334,
    'f1': 0.7368421052631579,
    'accuracy': None}