"""
Classifies: CHEBI:15693 aldose
"""
from rdkit import Chem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose (aldehydic parent sugars and their intramolecular hemiacetals).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldose, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of an aldehyde group or intramolecular hemiacetal
    aldehyde = Chem.MolFromSmarts('[CX3H1](=O)[#6]')
    hemiacetal = Chem.MolFromSmarts('C1(CO)OC(O)C1')
    if not mol.HasSubstructMatch(aldehyde) and not mol.HasSubstructMatch(hemiacetal):
        return False, "No aldehyde group or intramolecular hemiacetal found"

    # Check for the presence of at least two hydroxyl groups
    hydroxyl_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[CX4][OH]')))
    if hydroxyl_count < 2:
        return False, "Less than two hydroxyl groups found"

    return True, "Molecule is an aldose"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:15693',
                          'name': 'aldose',
                          'definition': 'Aldehydic parent sugars (polyhydroxy '
                                        'aldehydes H[CH(OH)]nC(=O)H, n >= 2) '
                                        'and their intramolecular hemiacetals.',
                          'parents': ['CHEBI:35381']},
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
    'num_true_positives': 5,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 18,
    'precision': 0.7142857142857143,
    'recall': 0.21739130434782608,
    'f1': 0.3333333333333333,
    'accuracy': None}