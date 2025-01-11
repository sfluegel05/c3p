"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion based on its SMILES string.
    A quaternary ammonium ion is defined as a derivative of ammonium (NH4(+)) where all four hydrogens
    are replaced by univalent (usually organyl) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quaternary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a nitrogen atom with four bonds and a positive charge
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen atom
            if atom.GetFormalCharge() == 1:  # Positive charge
                if len(atom.GetNeighbors()) == 4:  # Four bonds
                    # Check that all neighbors are carbon or other univalent groups
                    all_univalent = True
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() != 6:  # Not a carbon atom
                            all_univalent = False
                            break
                    if all_univalent:
                        return True, "Contains a nitrogen atom with four univalent groups and a positive charge"

    return False, "No quaternary ammonium ion pattern found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35267',
                          'name': 'quaternary ammonium ion',
                          'definition': 'A derivative of ammonium, NH4(+), in '
                                        'which all four of the hydrogens '
                                        'bonded to nitrogen have been replaced '
                                        'with univalent (usually organyl) '
                                        'groups.',
                          'parents': ['CHEBI:35267', 'CHEBI:35267']},
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