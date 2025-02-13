"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: nitrile (CHEBI:15718)
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile has a carbon atom triple-bonded to a nitrogen atom (RC≡N).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrile, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define nitrile SMARTS pattern (carbon triple-bonded to nitrogen)
    nitrile_pattern = Chem.MolFromSmarts("[CX2]#N")

    # Search for nitrile groups
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)
    num_nitrile_groups = len(nitrile_matches)

    if num_nitrile_groups == 0:
        return False, "No nitrile group (RC≡N) found"

    # Check that nitrile carbon is attached to at least one other carbon (R group)
    for match in nitrile_matches:
        nitrile_carbon_idx = match[0]
        nitrile_carbon = mol.GetAtomWithIdx(nitrile_carbon_idx)
        has_r_group = False
        for neighbor in nitrile_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != nitrile_carbon_idx:
                has_r_group = True
                break
        if not has_r_group:
            return False, "Nitrile carbon lacks R group (is hydrogen cyanide derivative)"

    return True, f"Contains {num_nitrile_groups} nitrile group(s) with R-C≡N structure"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:15718',
        'name': 'nitrile',
        'definition': 'A compound having the structure RC≡N; thus a C-substituted derivative of hydrocyanic acid, HC≡N. In systematic nomenclature, the suffix nitrile denotes the triply bound ≡N atom, not the carbon atom attached to it.',
        'parents': ['CHEBI:35342']
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