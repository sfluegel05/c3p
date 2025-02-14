"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
"""
Classifies: CHEBI:35733 11beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.
    An 11beta-hydroxy steroid is a steroid with a hydroxy group at position 11 in the beta configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for steroid backbone pattern
    # Allowing more flexibility in the pattern
    steroid_pattern = Chem.MolFromSmarts("[C@]12C[C@@H]3[C@H]([C@H]([C@]1(C)CCC4=CC(=O)CC[C@]34C)C)[C@@H]2CC")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Steroid backbone not found"

    # Look for 11beta-hydroxy group
    # Using a more flexible pattern
    hydroxy_pattern = Chem.MolFromSmarts("[C@H](O)[C@@]1([H])CCC2=CC(=O)CC[C@]12C")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "11beta-hydroxy group not found"

    # Check for typical molecular weight range (200-500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 500:
        return False, "Molecular weight outside typical range for 11beta-hydroxy steroids"

    return True, "Molecule contains steroid backbone with 11beta-hydroxy group"

# Example usage
print(is_11beta_hydroxy_steroid("CC(=O)C1CCC2C3CC(O)C4=CC(=O)CCC4(C)C3CCC12C"))  # True, '...'
print(is_11beta_hydroxy_steroid("C1CCCCC1"))  # False, '...'

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35733',
        'name': '11beta-hydroxy steroid',
        'definition': 'Any 11-hydroxy steroid in which the hydroxy group at position 11 has beta- configuration.',
        'parents': ['CHEBI:35619', 'CHEBI:35717']
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
    'num_true_positives': 141,
    'num_false_positives': 0,
    'num_true_negatives': 182423,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': 1.0
}