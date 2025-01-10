"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: CHEBI:38202 furanocoumarin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin is a furochromene that consists of a furan ring fused with a coumarin.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the coumarin pattern (benzene fused to pyrone)
    coumarin_pattern = Chem.MolFromSmarts("c1ccc2c(c1)C(=O)OC=C2")
    if not mol.HasSubstructMatch(coumarin_pattern):
        return False, "No coumarin structure found"

    # Define the furan pattern (five-membered ring with one oxygen)
    furan_pattern = Chem.MolFromSmarts("o1cccc1")
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found"

    # Check if the furan ring is fused to the coumarin structure
    # We look for a pattern where the furan ring is fused to the coumarin
    fused_pattern = Chem.MolFromSmarts("c1ccc2c(c1)C(=O)OC=C2.o1cccc1")
    if not mol.HasSubstructMatch(fused_pattern):
        return False, "Furan ring is not fused to the coumarin structure"

    return True, "Contains a furan ring fused to a coumarin structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38202',
                          'name': 'furanocoumarin',
                          'definition': 'Any furochromene that consists of a '
                                        'furan ring fused with a coumarin. The '
                                        'fusion may occur in different ways in '
                                        'give several isomers.'},
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