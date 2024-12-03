"""
Classifies: CHEBI:33296 alkali metal molecular entity
"""
from rdkit import Chem

def is_alkali_metal_molecular_entity(smiles: str):
    """
    Determines if a molecule is an alkali metal molecular entity (contains one or more atoms of an alkali metal).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkali metal molecular entity, False otherwise
        str: Reason for classification
    """
    alkali_metals = {'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr'}
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    for atom in mol.GetAtoms():
        if atom.GetSymbol() in alkali_metals:
            return True, f"Contains alkali metal atom: {atom.GetSymbol()}"

    return False, "No alkali metal atoms found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33296',
                          'name': 'alkali metal molecular entity',
                          'definition': 'A molecular entity containing one or '
                                        'more atoms of an alkali metal.',
                          'parents': ['CHEBI:33674']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 29-30: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}