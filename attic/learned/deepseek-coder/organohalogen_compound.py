"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: CHEBI:37141 organohalogen compound
"""
from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.
    An organohalogen compound contains at least one carbon-halogen bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organohalogen compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define halogens
    halogens = {'F', 'Cl', 'Br', 'I', 'At'}

    # Iterate through atoms to find halogens bonded to carbon
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in halogens:
            # Check if the halogen is bonded to at least one carbon atom
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    return True, f"Contains a carbon-halogen bond with {atom.GetSymbol()}"

    # If no carbon-halogen bond is found
    return False, "No carbon-halogen bond found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37141',
                          'name': 'organohalogen compound',
                          'definition': 'A compound containing at least one '
                                        'carbon-halogen bond (where X is a '
                                        'halogen atom).',
                          'parents': ['CHEBI:24431', 'CHEBI:37141']},
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