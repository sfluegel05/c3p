"""
Classifies: CHEBI:17792 organohalogen compound
"""
from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.
    An organohalogen compound contains at least one carbon-halogen bond (where halogen is F, Cl, Br, or I).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an organohalogen compound, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for carbon-halogen bond
    halogen_pattern = Chem.MolFromSmarts("[#6]-[F,Cl,Br,I]")
    if halogen_pattern is None:
        return False, "Invalid SMARTS pattern for halogens"
    
    # Check if the molecule contains at least one carbon-halogen bond
    if mol.HasSubstructMatch(halogen_pattern):
        return True, "Contains at least one carbon-halogen bond"
    else:
        return False, "No carbon-halogen bond found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:24432',
        'name': 'organohalogen compound',
        'definition': 'A compound containing at least one carbon-halogen bond (where X is a halogen atom).',
        'parents': []
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
    'attempt': 0,
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