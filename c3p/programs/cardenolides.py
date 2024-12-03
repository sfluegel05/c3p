"""
Classifies: CHEBI:74634 cardenolides
"""
from rdkit import Chem

def is_cardenolides(smiles: str):
    """
    Determines if a molecule is a cardenolide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardenolide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for C23 steroid backbone
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("C1CCC2C1(CCC3C2CCC4C3(CCC5C4C(CC(C5)C=O)O)C)C")):
        return False, "Molecule does not have a C23 steroid backbone"

    # Check for five-membered lactone ring at C-17
    lactone_ring = Chem.MolFromSmarts("C1=CCOC1=O")
    if not mol.HasSubstructMatch(lactone_ring):
        return False, "Molecule does not have a five-membered lactone ring"

    return True, "Molecule is a cardenolide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:74634',
                          'name': 'cardenolides',
                          'definition': 'Any steroid lactone that is a  C23 '
                                        'steroid with a five-membered lactone '
                                        'ring at C-17 and its substituted '
                                        'derivatives. They form the aglycone '
                                        'constituents of cardiac glycosides.',
                          'parents': ['CHEBI:26766', 'CHEBI:50523']},
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
    'num_true_negatives': 11,
    'num_false_negatives': 11,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}