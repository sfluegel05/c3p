"""
Classifies: CHEBI:24866 salt
"""
from rdkit import Chem

def is_salt(smiles: str):
    """
    Determines if a molecule is a salt (an assembly of cations and anions).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a salt, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of cations and anions
    cations = [atom for atom in mol.GetAtoms() if atom.GetFormalCharge() > 0]
    anions = [atom for atom in mol.GetAtoms() if atom.GetFormalCharge() < 0]

    if cations and anions:
        return True, "Molecule contains both cations and anions"
    else:
        return False, "Molecule does not contain both cations and anions"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24866',
                          'name': 'salt',
                          'definition': 'A salt is an assembly of cations and '
                                        'anions.',
                          'parents': ['CHEBI:37577']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 57-58: malformed \\N character escape (<string>, line '
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