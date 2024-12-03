"""
Classifies: CHEBI:50753 isoflavonoid
"""
from rdkit import Chem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 1-benzopyran core
    benzopyran_smarts = "c1cc2occc2cc1"
    benzopyran = Chem.MolFromSmarts(benzopyran_smarts)
    if not mol.HasSubstructMatch(benzopyran):
        return False, "No 1-benzopyran core found"

    # Check for aryl substituent at position 3
    aryl_substituent_smarts = "c1cc2occc2cc1-c3ccccc3"
    aryl_substituent = Chem.MolFromSmarts(aryl_substituent_smarts)
    if not mol.HasSubstructMatch(aryl_substituent):
        return False, "No aryl substituent at position 3 found"

    return True, "Molecule is an isoflavonoid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50753',
                          'name': 'isoflavonoid',
                          'definition': 'Any 1-benzopyran with an aryl '
                                        'substituent at position 3. The term '
                                        'was originally restricted to natural '
                                        'products, but is now also used to '
                                        'describe semi-synthetic and fully '
                                        'synthetic compounds.',
                          'parents': [   'CHEBI:26004',
                                         'CHEBI:38443',
                                         'CHEBI:72544']},
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
    'num_false_negatives': 63,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}