"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: CHEBI:138100
beta-D-glucoside
Definition: Any D-glucoside in which the anomeric centre has beta-configuration.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside is a glycoside in which the glucose moiety is in the D-configuration
    and the anomeric carbon (C1) has the beta-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the beta-D-glucoside SMARTS pattern
    # This pattern matches beta-D-glucopyranose attached via the anomeric oxygen
    beta_D_glucoside_smarts = """
    [C@@H]1([O][#6])[O][C@@H]([C@@H]([C@H]([C@H]1O)O)O)O
    """

    # Remove whitespace and newlines from SMARTS pattern
    beta_D_glucoside_smarts = beta_D_glucoside_smarts.replace('\n', '').replace(' ', '')

    pattern = Chem.MolFromSmarts(beta_D_glucoside_smarts)
    if pattern is None:
        return False, "Invalid SMARTS pattern for beta-D-glucoside"

    # Check for substructure match
    if mol.HasSubstructMatch(pattern):
        return True, "Contains beta-D-glucoside moiety with beta-configuration at the anomeric carbon"
    else:
        return False, "Does not contain beta-D-glucoside moiety with correct stereochemistry"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:138100',
        'name': 'beta-D-glucoside',
        'definition': 'Any D-glucoside in which the anomeric centre has beta-configuration.',
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
    'attempt': 1,
    'success': None,
    'best': None,
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