"""
Classifies: CHEBI:33521 metal atom
"""
"""
Classifies: CHEBI:33521 metal atom
"""
from rdkit import Chem

def is_metal_atom(smiles: str):
    """
    Determines if a molecule is a metal atom based on its SMILES string.
    A metal atom is defined as an atom of an element that exhibits typical metallic properties.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a metal atom, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule consists of a single atom
    if mol.GetNumAtoms() != 1:
        return False, "Molecule must consist of a single atom"

    # Get the atomic number of the atom
    atom = mol.GetAtomWithIdx(0)
    atomic_number = atom.GetAtomicNum()

    # List of atomic numbers corresponding to metals
    metal_atomic_numbers = set(range(3, 4)) | set(range(11, 13)) | set(range(19, 32)) | \
                           set(range(37, 51)) | set(range(55, 84)) | set(range(87, 104)) | \
                           set(range(104, 119))

    # Check if the atomic number is in the list of metals
    if atomic_number in metal_atomic_numbers:
        return True, f"Atom with atomic number {atomic_number} is a metal"
    else:
        return False, f"Atom with atomic number {atomic_number} is not a metal"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33521',
                          'name': 'metal atom',
                          'definition': 'An atom of an element that exhibits '
                                        'typical metallic properties, being '
                                        'typically shiny, with high electrical '
                                        'and thermal conductivity.',
                          'parents': ['CHEBI:24835', 'CHEBI:24836']},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}