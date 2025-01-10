"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
"""
Classifies: CHEBI:37653 3-oxo-Delta(4) steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    A 3-oxo-Delta(4) steroid has a steroid backbone with a 3-oxo group and a C=C double bond at the alpha,beta position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(4) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid backbone (fused ring system with at least 17 carbons)
    steroid_pattern = Chem.MolFromSmarts("[C]1[C]2[C]3[C]4[C]5[C]6[C]7[C]8[C]9[C]%10[C]%11[C]%12[C]%13[C]%14[C]%15[C]%16[C]%17")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for 3-oxo group (carbonyl at position 3)
    oxo_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if not any(match[0].GetAtomWithIdx(0).GetIdx() == 3 for match in oxo_matches):
        return False, "No 3-oxo group found"

    # Check for C=C double bond at position 4 (alpha,beta to the 3-oxo group)
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if not any(match[0].GetAtomWithIdx(0).GetIdx() == 4 for match in double_bond_matches):
        return False, "No C=C double bond at position 4 found"

    return True, "Contains steroid backbone with 3-oxo group and C=C double bond at position 4"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37653',
                          'name': '3-oxo-Delta(4) steroid',
                          'definition': 'A 3-oxo steroid conjugated to a C=C double bond at the alpha,beta position.',
                          'parents': ['CHEBI:37653', 'CHEBI:37653']},
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
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}