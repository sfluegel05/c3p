"""
Classifies: CHEBI:36885 20-oxo steroid
"""
from rdkit import Chem

def is_20_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 20-oxo steroid, defined as 'An oxo steroid carrying an oxo group at position 20.'.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 20-oxo steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a steroid backbone
    steroid_smarts = "[#6]1[#6][#6][#6]2[#6][#6][#6]3[#6]([#6][#6][#6]4[#6](=[#8])[#6][#6][#6][#6]4[#6]3[#6][#6]2[#6]1)"
    steroid = Chem.MolFromSmarts(steroid_smarts)
    
    if not mol.HasSubstructMatch(steroid):
        return False, "Molecule does not have a steroid backbone"

    # Define the SMARTS pattern for a 20-oxo group
    oxo_20_smarts = "[#6]1[#6][#6][#6]2[#6][#6][#6]3[#6]([#6][#6][#6]4[#6](=[#8])[#6][#6][#6][#6]4[#6]3[#6][#6]2[#6]1)C(=O)"
    oxo_20 = Chem.MolFromSmarts(oxo_20_smarts)
    
    if mol.HasSubstructMatch(oxo_20):
        return True, "Molecule is a 20-oxo steroid"
    else:
        return False, "Molecule does not have an oxo group at position 20"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36885',
                          'name': '20-oxo steroid',
                          'definition': 'An oxo steroid carrying an oxo group '
                                        'at position 20.',
                          'parents': ['CHEBI:35789']},
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
    'num_true_negatives': 19,
    'num_false_negatives': 19,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}