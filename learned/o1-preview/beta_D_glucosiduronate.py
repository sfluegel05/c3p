"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: CHEBI:59687 beta-D-glucosiduronate
"""
from rdkit import Chem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronate based on its SMILES string.
    A beta-D-glucosiduronate is a compound containing a beta-D-glucuronic acid moiety
    attached via a glycosidic bond, with a deprotonated carboxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for beta-D-glucuronic acid moiety connected via beta-glycosidic bond
    # Anomeric carbon (C1) has beta configuration ([C@H]) and is connected via oxygen to any carbon ([O][#6])
    # The ring is a pyranose with hydroxyl groups at positions 2,3,4
    # The C6 position is oxidized to a carboxylate (C(=O)[O-] or C(=O)O)
    beta_D_glucuronide_smarts = '[C@H]1([O][#6])O[C@@H]([C@@H](O)[C@@H](O)[C@H](O1)C(=O)[O-,O])'

    pattern = Chem.MolFromSmarts(beta_D_glucuronide_smarts)
    if pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check for substructure match
    if not mol.HasSubstructMatch(pattern):
        return False, "Beta-D-glucuronic acid moiety not found with correct stereochemistry"

    return True, "Contains beta-D-glucuronic acid moiety connected via beta-glycosidic bond"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:59687',
        'name': 'beta-D-glucosiduronate',
        'definition': 'A carbohydrate acid derivative anion obtained by deprotonation of the carboxy group of any beta-D-glucosiduronic acid; major species at pH 7.3.',
        'parents': ['CHEBI:17144', 'CHEBI:58395']
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
}