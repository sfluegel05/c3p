"""
Classifies: CHEBI:36683 organochlorine compound
"""
from rdkit import Chem

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound (contains at least one carbon-chlorine bond).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organochlorine compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through all the atoms in the molecule
    for atom in mol.GetAtoms():
        # Check if the atom is a chlorine
        if atom.GetSymbol() == 'Cl':
            # Check if the chlorine is bonded to a carbon
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    return True, "Contains at least one carbon-chlorine bond"
    
    return False, "No carbon-chlorine bonds found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36683',
                          'name': 'organochlorine compound',
                          'definition': 'An organochlorine compound is a '
                                        'compound containing at least one '
                                        'carbon-chlorine bond.',
                          'parents': ['CHEBI:17792', 'CHEBI:23117']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 10-11: malformed \\N character escape (<string>, line '
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