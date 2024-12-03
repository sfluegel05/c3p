"""
Classifies: CHEBI:37142 organoiodine compound
"""
from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound (contains at least one carbon-iodine bond).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organoiodine compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one carbon-iodine bond
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if (atom1.GetSymbol() == 'C' and atom2.GetSymbol() == 'I') or (atom1.GetSymbol() == 'I' and atom2.GetSymbol() == 'C'):
            return True, "Contains at least one carbon-iodine bond"

    return False, "No carbon-iodine bonds found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37142',
                          'name': 'organoiodine compound',
                          'definition': 'An organoiodine compound is a '
                                        'compound containing at least one '
                                        'carbon-iodine bond.',
                          'parents': ['CHEBI:17792', 'CHEBI:24860']},
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
    'num_true_positives': 11,
    'num_false_positives': 0,
    'num_true_negatives': 11,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}