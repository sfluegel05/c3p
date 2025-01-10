"""
Classifies: CHEBI:72544 flavonoids
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids are characterized by a C6-C3-C6 carbon skeleton with a benzene ring
    attached to a heterocyclic benzopyran ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for flavonoid benzopyran backbone pattern, C6-C3-C6 structure 
    flavonoid_pattern = Chem.MolFromSmarts("c1cc(O)c2c(c1)-c3c(ccc(o3))-[C6]")
    if not mol.HasSubstructMatch(flavonoid_pattern):
        return False, "No flavonoid backbone found"

    return True, "Contains flavonoid backbone (C6-C3-C6 structure with benzopyran)"

__metadata__ = {   'chemical_class': {   'id': 'None',
                          'name': 'flavonoids',
                          'definition': 'Compounds characterized by a C6-C3-C6 carbon skeleton with a benzene ring '
                                        'attached to a heterocyclic benzopyran ring.',
                          'parents': ['CHEBI:24013']},
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
    'success': None,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'num_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': 0.0}