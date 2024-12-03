"""
Classifies: CHEBI:33318 main group element atom
"""
from rdkit import Chem

def is_main_group_element_atom(smiles: str):
    """
    Determines if a molecule is a main group element atom (found in the s- and p- blocks of the periodic table).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a main group element atom, False otherwise
        str: Reason for classification
    """
    main_group_elements = {
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 
        'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 
        'K', 'Ca', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 
        'Rb', 'Sr', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 
        'Cs', 'Ba', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 
        'Fr', 'Ra'
    }

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    atoms = mol.GetAtoms()
    if len(atoms) != 1:
        return False, "Not a single atom"

    atom = atoms[0]
    symbol = atom.GetSymbol()

    if symbol in main_group_elements:
        return True, f"Atom is a main group element: {symbol}"
    else:
        return False, f"Atom is not a main group element: {symbol}"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33318',
                          'name': 'main group element atom',
                          'definition': 'An atom belonging to one of the main '
                                        'groups (found in the s- and p- '
                                        'blocks) of the periodic table.',
                          'parents': ['CHEBI:33250']},
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
    'num_true_positives': 14,
    'num_false_positives': 0,
    'num_true_negatives': 14,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}