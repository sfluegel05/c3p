"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
"""
Classifies: 4'-hydroxyflavanones
Definition: Any hydroxyflavanone having a hydroxy substituent located at position 4'.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    A 4'-hydroxyflavanone has a flavanone backbone with a hydroxy group at the 4' position of the phenyl ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 4'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavanone backbone pattern (2,3-dihydro-2-phenylchromen-4-one)
    flavanone_pattern = Chem.MolFromSmarts("[O;H0]c1ccc2C(=O)CC(Oc2c1)c1ccc([OH])cc1")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone backbone found"

    # Check for the hydroxy group at the 4' position (4th carbon of the phenyl ring)
    hydroxy_4_prime_pattern = Chem.MolFromSmarts("[OH]c1ccc(cc1)C2CC(=O)c3ccccc3O2")
    if not mol.HasSubstructMatch(hydroxy_4_prime_pattern):
        return False, "No hydroxy group at the 4' position"

    return True, "Contains flavanone backbone with a hydroxy group at the 4' position"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:XXXXX',
                          'name': '4\'-hydroxyflavanones',
                          'definition': 'Any hydroxyflavanone having a hydroxy substituent located at position 4\'.',
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': 0.0}