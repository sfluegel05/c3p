"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: beta-D-glucoside
"""
from rdkit import Chem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside is a molecule containing a beta-D-glucose unit where the anomeric carbon (C1) is connected via an O-glycosidic bond in the beta configuration.

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

    # Define the SMARTS pattern for beta-D-glucoside
    beta_D_glucoside_smarts = """
    [C@H]1([O][#6])[O][C@@H]([C@@H]([C@H]([C@H]1O)O)O)O
    """
    pattern = Chem.MolFromSmarts(beta_D_glucoside_smarts.strip())
    if pattern is None:
        return False, "Failed to parse SMARTS pattern"

    # Check for beta-D-glucoside substructure
    if mol.HasSubstructMatch(pattern):
        return True, "Contains beta-D-glucoside moiety"
    else:
        return False, "Does not contain beta-D-glucoside moiety"

__metadata__ = {   'chemical_class': {   'id': None,
                              'name': 'beta-D-glucoside',
                              'definition': 'Any D-glucoside in which the anomeric centre has beta-configuration.',
                              'parents': []},
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
        'accuracy': None}