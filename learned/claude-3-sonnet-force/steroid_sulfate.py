"""
Classifies: CHEBI:16158 steroid sulfate
"""
"""
Classifies: CHEBI:38866 steroid sulfate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    A steroid sulfate is a sulfuric ester obtained by the formal condensation of 
    a hydroxy group of any steroid with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid sulfate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for steroid backbone pattern (tetracyclic ring system with specific ring junctions)
    steroid_pattern = Chem.MolFromSmarts(
        "[C@]12[C@@]([C@](C[C@@]1[C@]([C@@]([C@]3([C@]([C@]([C@@]2(C)C)([H])CC3)([H])CC)([H])C)([H])[H])CC)([H])[H]")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Look for sulfate group (-O-S(=O)(=O)-O)
    sulfate_pattern = Chem.MolFromSmarts("[OX2][SX4](=[OX1])(=[OX1])[OX2]")
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)
    if not sulfate_matches:
        return False, "No sulfate group found"

    # Check for common functional groups found in steroid sulfates
    functional_groups = ["[C@H](CC)", "[C@H](O)", "[C@H](CCC)", "[C@H](CC=O)"]
    functional_group_present = any(mol.HasSubstructMatch(Chem.MolFromSmarts(fg)) for fg in functional_groups)
    if not functional_group_present:
        return False, "Missing common functional groups found in steroid sulfates"

    # Count rotatable bonds to verify steroid-like rigidity
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable > 10:
        return False, "Too many rotatable bonds for a steroid-like structure"

    # Check molecular weight - steroid sulfates typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for steroid sulfate"

    return True, "Contains steroid backbone and sulfate group, with common functional groups found in steroid sulfates"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:38866',
        'name': 'steroid sulfate',
        'definition': 'A sulfuric ester obtained by the formal condensation of a hydroxy group of any steroid with sulfuric acid.',
        'parents': ['CHEBI:35695', 'CHEBI:24347']
    },
    'config': {
        'llm_model_name': 'claude-sonnet',
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
    'success': False,
    'best': False,
    'error': '',
    'stdout': None,
    'num_true_positives': 85,
    'num_false_positives': 0,
    'num_true_negatives': 182419,
    'num_false_negatives': 80,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.515151515151515,
    'f1': 0.6802721088435374,
    'accuracy': 0.9995632778394829
}