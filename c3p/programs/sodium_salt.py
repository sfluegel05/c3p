"""
Classifies: CHEBI:26714 sodium salt
"""
from rdkit import Chem

def is_sodium_salt(smiles: str):
    """
    Determines if a molecule is a sodium salt (any alkali metal salt having sodium(1+) as the cation).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sodium salt, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    sodium_found = False
    anion_found = False

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "Na" and atom.GetFormalCharge() == 1:
            sodium_found = True
        elif atom.GetFormalCharge() == -1:
            anion_found = True

    if sodium_found and anion_found:
        return True, "Molecule is a sodium salt"
    elif not sodium_found:
        return False, "No sodium cation found"
    elif not anion_found:
        return False, "No anion found"
    
    return None, None


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26714',
                          'name': 'sodium salt',
                          'definition': 'Any alkali metal salt having '
                                        'sodium(1+) as the cation.',
                          'parents': ['CHEBI:26712', 'CHEBI:35479']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 23-24: malformed \\N character escape (<string>, line '
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