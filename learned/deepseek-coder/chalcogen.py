"""
Classifies: CHEBI:33303 chalcogen
"""
"""
Classifies: CHEBI:33302 chalcogen
"""
from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a molecule is a chalcogen based on its SMILES string.
    A chalcogen is any p-block element belonging to the group 16 family of the periodic table.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcogen, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is a single atom
    if mol.GetNumAtoms() != 1:
        return False, "Not a single atom"

    # Get the atomic number of the atom
    atom = mol.GetAtomWithIdx(0)
    atomic_num = atom.GetAtomicNum()

    # List of atomic numbers for chalcogens: O (8), S (16), Se (34), Te (52), Po (84)
    chalcogen_atomic_nums = {8, 16, 34, 52, 84}

    if atomic_num in chalcogen_atomic_nums:
        # Ensure the atom is neutral (no charge)
        if atom.GetFormalCharge() == 0:
            return True, f"Element with atomic number {atomic_num} is a chalcogen"
        else:
            return False, f"Element with atomic number {atomic_num} has a charge and is not a neutral chalcogen atom"
    else:
        return False, f"Element with atomic number {atomic_num} is not a chalcogen"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33302',
                          'name': 'chalcogen',
                          'definition': 'Any p-block element belonging to the '
                                        'group 16 family of the periodic '
                                        'table.',
                          'parents': ['CHEBI:33302']},
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