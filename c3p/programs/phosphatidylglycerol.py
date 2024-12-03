"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
from rdkit import Chem

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylglycerol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMILES pattern for a phosphatidyl group
    phosphatidyl_pattern = Chem.MolFromSmarts('O=P(O)(O)OCC(CO)O')
    if phosphatidyl_pattern is None:
        return False, "Invalid phosphatidyl pattern"

    # Define the SMILES pattern for glycerol
    glycerol_pattern = Chem.MolFromSmarts('OCC(CO)O')
    if glycerol_pattern is None:
        return False, "Invalid glycerol pattern"

    # Check if the molecule contains the phosphatidyl group
    if not mol.HasSubstructMatch(phosphatidyl_pattern):
        return False, "No phosphatidyl group found"

    # Check if the molecule contains the glycerol group
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol group found"

    # Check if the phosphatidyl group is attached to the glycerol group
    phosphatidyl_matches = mol.GetSubstructMatches(phosphatidyl_pattern)
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    for phosphatidyl_match in phosphatidyl_matches:
        for glycerol_match in glycerol_matches:
            if set(phosphatidyl_match).intersection(set(glycerol_match)):
                return True, "Molecule is a phosphatidylglycerol"

    return False, "Phosphatidyl group is not attached to glycerol group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17517',
                          'name': 'phosphatidylglycerol',
                          'definition': 'A glycerophosphoglycerol that is '
                                        'glycerol in which the hydrogen of one '
                                        'of the primary hydroxy groups has '
                                        'been replaced by a phosphatidyl '
                                        'group.',
                          'parents': ['CHEBI:24360']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 50,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 0,
    'precision': 0.9803921568627451,
    'recall': 1.0,
    'f1': 0.99009900990099,
    'accuracy': None}