"""
Classifies: CHEBI:46633 carbapenems
"""
from rdkit import Chem

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    A carbapenem is a beta-lactam antibiotic with a carbapenem skeleton,
    which is a fused beta-lactam ring and a five-membered ring with an exocyclic double bond,
    variously substituted at positions 3, 4, and 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbapenem, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Carbapenem core SMARTS pattern
    # This pattern captures:
    # - Fused bicyclic system: four-membered beta-lactam ring fused to a five-membered ring
    # - Beta-lactam ring: four-membered ring with N and carbonyl group
    # - Five-membered ring: with exocyclic double bond
    # - Substituents allowed at positions 3, 4, and 6

    carbapenem_smarts = '''
    [C;R1]=C1[C;R1][C;R1][C;R1]2[N;R1][C;R1](=O)[C;R1]12
    '''
    carbapenem_pattern = Chem.MolFromSmarts(carbapenem_smarts)
    if carbapenem_pattern is None:
        return False, "Error in defining carbapenem SMARTS pattern"

    # Check if molecule contains the carbapenem core
    if not mol.HasSubstructMatch(carbapenem_pattern):
        return False, "Does not contain carbapenem core structure"

    # Check for exocyclic double bond on five-membered ring
    # The exocyclic double bond is between a ring carbon and an external carbon
    exocyclic_double_bond_pattern = Chem.MolFromSmarts('[C;R][C;R]=[C;R0]')
    if not mol.HasSubstructMatch(exocyclic_double_bond_pattern):
        return False, "Missing exocyclic double bond on five-membered ring"

    return True, "Contains carbapenem core structure"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:23066',
        'name': 'carbapenems',
        'definition': 'The class of beta-lactam antibiotics whose members have a carbapenem skeleton which is variously substituted at positions 3, 4, and 6.',
        'parents': ['CHEBI:19258', 'CHEBI:26836']
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}