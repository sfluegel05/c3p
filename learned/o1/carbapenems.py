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

    # Define the carbapenem core structure using SMARTS
    # This pattern captures:
    # - A fused bicyclic ring system with a four-membered beta-lactam ring and a five-membered ring
    # - Beta-lactam ring contains nitrogen and carbonyl group
    # - Five-membered ring with exocyclic double bond

    carbapenem_smarts = '''
    [$([C;R3]-1-[C;R3]-[N;R3]-[C;R3](=O)-1)]
    [$([C;R5]-2=[C;R5]-[C;R5]-[C;R5]-[C;R5]-2)]
    '''

    # Remove line breaks and whitespace
    carbapenem_smarts = ''.join(carbapenem_smarts.split())

    # Combine the two patterns to represent the fused ring system
    # The exocyclic double bond is attached to the five-membered ring
    carbapenem_pattern = Chem.MolFromSmarts(carbapenem_smarts)
    if carbapenem_pattern is None:
        return False, "Error in defining carbapenem SMARTS pattern"

    # Check if molecule contains the carbapenem core
    if not mol.HasSubstructMatch(carbapenem_pattern):
        return False, "Does not contain carbapenem core structure"

    # Additional check for exocyclic double bond at the five-membered ring
    # Look for C=C double bond attached to the five-membered ring
    exocyclic_double_bond = Chem.MolFromSmarts('[C;R5]=[C;R0]')
    if not mol.HasSubstructMatch(exocyclic_double_bond):
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
    'attempt': 1,
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