"""
Classifies: CHEBI:50994 primary amino compound
"""
from rdkit import Chem

def is_primary_amino_compound(smiles: str):
    """
    Determines if a molecule is a primary amino compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amino compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through all atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Check if the atom is nitrogen
            num_H = sum([1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 1])
            num_C = sum([1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6])
            if num_H == 2 and num_C == 1:  # Check for NH2 group attached to one carbon
                return True, "Primary amino group (NH2) attached to one carbon found"
    
    return False, "No primary amino group (NH2) attached to one carbon found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50994',
                          'name': 'primary amino compound',
                          'definition': 'A compound formally derived from '
                                        'ammonia by replacing one hydrogen '
                                        'atom by an organyl group.',
                          'parents': ['CHEBI:50047']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 11-12: malformed \\N character escape (<string>, line '
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