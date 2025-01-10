"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
"""
Classifies: CHEBI:XXXXX 3'-hydroxyflavanones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    A 3'-hydroxyflavanone is a flavanone with a hydroxyl group at the 3' position of the phenyl ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavanone backbone pattern
    flavanone_pattern = Chem.MolFromSmarts("[C@@H]1CC(=O)c2ccccc2O1")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone backbone found"

    # Define the 3'-hydroxyflavanone pattern
    # The pattern matches the flavanone backbone with a hydroxyl group at the 3' position of the phenyl ring
    hydroxyflavanone_pattern = Chem.MolFromSmarts("[C@@H]1CC(=O)c2cc(O)ccc2O1")
    if not mol.HasSubstructMatch(hydroxyflavanone_pattern):
        return False, "No hydroxyl group at the 3' position of the phenyl ring"

    return True, "Contains flavanone backbone with a hydroxyl group at the 3' position of the phenyl ring"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:XXXXX',
                          'name': '3\'-hydroxyflavanones',
                          'definition': 'Any hydroxyflavanone with a hydroxy '
                                        'substituent at position 3\' of the '
                                        'phenyl ring.',
                          'parents': ['CHEBI:XXXXX', 'CHEBI:XXXXX']},
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