"""
Classifies: CHEBI:23666 diamine
"""
from rdkit import Chem

def is_diamine(smiles: str):
    """
    Determines if a molecule is a diamine (contains two amino groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all amino groups (NH2, NH, or N)
    amino_groups = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'N']

    if len(amino_groups) == 2:
        return True, "Contains exactly two amino groups"
    elif len(amino_groups) < 2:
        return False, "Contains less than two amino groups"
    else:
        return False, "Contains more than two amino groups"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23666',
                          'name': 'diamine',
                          'definition': 'Any polyamine that contains two amino '
                                        'groups.',
                          'parents': ['CHEBI:88061']},
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
    'num_false_positives': 0,
    'num_true_negatives': 13,
    'num_false_negatives': 8,
    'precision': 1.0,
    'recall': 0.38461538461538464,
    'f1': 0.5555555555555556,
    'accuracy': None}