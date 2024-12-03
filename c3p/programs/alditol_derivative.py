"""
Classifies: CHEBI:63423 alditol derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alditol_derivative(smiles: str):
    """
    Determines if a molecule is an alditol derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a glycerol backbone or similar structures
    glycerol_smarts = Chem.MolFromSmarts("C(CO)CO")
    if not mol.HasSubstructMatch(glycerol_smarts):
        return False, "No glycerol backbone found"

    # Check for phosphate groups
    phosphate_smarts = Chem.MolFromSmarts("P(=O)(O)(O)")
    if mol.HasSubstructMatch(phosphate_smarts):
        return True, "Contains glycerol backbone and phosphate group"

    # Check for acyl groups
    acyl_smarts = Chem.MolFromSmarts("C(=O)O")
    if mol.HasSubstructMatch(acyl_smarts):
        return True, "Contains glycerol backbone and acyl group"

    return False, "Does not match known alditol derivative structures"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:63423',
                          'name': 'alditol derivative',
                          'definition': 'A carbohydrate derivative that is '
                                        'formally obtained from an alditol.',
                          'parents': ['CHEBI:63299']},
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
    'num_true_positives': 85,
    'num_false_positives': 5,
    'num_true_negatives': 15,
    'num_false_negatives': 5,
    'precision': 0.9444444444444444,
    'recall': 0.9444444444444444,
    'f1': 0.9444444444444444,
    'accuracy': None}