"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a phosphate group
    phosphate = Chem.MolFromSmarts("P(=O)(O)(O)O")
    if not mol.HasSubstructMatch(phosphate):
        return False, "No phosphate group found"

    # Check for the presence of a prenol backbone
    prenol = Chem.MolFromSmarts("C=C(C)CC")
    if not mol.HasSubstructMatch(prenol):
        return False, "No prenol backbone found"

    # Check for the formal condensation of the terminal allylic hydroxy group with phosphoric acid
    prenol_phosphate = Chem.MolFromSmarts("C=C(C)CCOP(=O)(O)O")
    if not mol.HasSubstructMatch(prenol_phosphate):
        return False, "No formal condensation of the terminal allylic hydroxy group with phosphoric acid found"

    return True, "Molecule is a polyprenol phosphate"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16460',
                          'name': 'polyprenol phosphate',
                          'definition': 'A prenol phosphate resulting from the '
                                        'formal condensation of the terminal '
                                        'allylic hydroxy group of a polyprenol '
                                        'with 1 mol eq. of phosphoric acid.',
                          'parents': ['CHEBI:26250', 'CHEBI:26875']},
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