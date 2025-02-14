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
    steroid_pattern1 = Chem.MolFromSmarts("[C@]12[C@@]([C@](C[C@@]1[C@]([C@@]([C@]3([C@]([C@]([C@@]2(C)C)([H])CC3)([H])CC)([H])[H])CC)([H])[H]")
    steroid_pattern2 = Chem.MolFromSmarts("[C@]12[C@@]([C@](C[C@@]1[C@]([C@@]([C@]3([C@]([C@]([C@@]2(C)C)([H])C(C)C3)([H])C)([H])[H])CC)([H])[H]")
    if not (mol.HasSubstructMatch(steroid_pattern1) or mol.HasSubstructMatch(steroid_pattern2)):
        return False, "No steroid backbone found"

    # Look for sulfate group (-O-S(=O)(=O)-O)
    sulfate_pattern = Chem.MolFromSmarts("[OX2][SX4](=[OX1])(=[OX1])[OX2]")
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)
    if not sulfate_matches:
        return False, "No sulfate group found"

    # Check molecular weight - steroid sulfates typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for steroid sulfate"

    return True, "Contains steroid backbone and sulfate group"


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
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 154,
    'num_false_positives': 0,
    'num_true_negatives': 182419,
    'num_false_negatives': 11,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.9333333333333333,
    'f1': 0.9655172413793103,
    'accuracy': 0.9994009567018254
}