"""
Classifies: CHEBI:29089 1,2-diacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem


def is_1_2_diacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a 1,2-diacyl-sn-glycerol 3-phosphate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1,2-diacyl-sn-glycerol 3-phosphate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the glycerol backbone
    glycerol = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol):
        return False, "No glycerol backbone found"

    # Check for the presence of the phosphate group
    phosphate = Chem.MolFromSmarts("COP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate):
        return False, "No phosphate group found"

    # Check for acyl groups at the 1- and 2- positions
    acyl1 = Chem.MolFromSmarts("OCC(OC(=O)C)CO")
    acyl2 = Chem.MolFromSmarts("OCC(O)COC(=O)C")
    if not (mol.HasSubstructMatch(acyl1) and mol.HasSubstructMatch(acyl2)):
        return False, "No acyl groups at the 1- and 2- positions"

    return True, "Molecule is a 1,2-diacyl-sn-glycerol 3-phosphate"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:29089',
                          'name': '1,2-diacyl-sn-glycerol 3-phosphate',
                          'definition': 'An sn-glycerol 3-phosphate compound '
                                        'having unspecified O-acyl groups at '
                                        'the 1- and 2-positions.',
                          'parents': ['CHEBI:16337', 'CHEBI:26706']},
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
    'num_true_positives': 43,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 2,
    'precision': 0.9772727272727273,
    'recall': 0.9555555555555556,
    'f1': 0.9662921348314608,
    'accuracy': None}