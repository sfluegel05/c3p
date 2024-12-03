"""
Classifies: CHEBI:16337 phosphatidic acid
"""
from rdkit import Chem

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the required substructures for phosphatidic acid
    glycerol_substructure = Chem.MolFromSmarts("OCC(O)CO")
    phosphate_substructure = Chem.MolFromSmarts("COP(=O)(O)O")
    ester_substructure = Chem.MolFromSmarts("C(=O)O")

    # Check for the presence of glycerol backbone
    if not mol.HasSubstructMatch(glycerol_substructure):
        return False, "No glycerol backbone found"

    # Check for the presence of phosphate group
    if not mol.HasSubstructMatch(phosphate_substructure):
        return False, "No phosphate group found"

    # Check for the presence of two ester bonds
    ester_matches = mol.GetSubstructMatches(ester_substructure)
    if len(ester_matches) < 2:
        return False, "Less than two ester bonds found"

    return True, "Phosphatidic acid detected"



__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16337',
                          'name': 'phosphatidic acid',
                          'definition': 'A derivative of glycerol in which one '
                                        'hydroxy group, commonly but not '
                                        'necessarily primary, is esterified '
                                        'with phosphoric acid and the other '
                                        'two are esterified with fatty acids.',
                          'parents': ['CHEBI:37739']},
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
    'num_true_positives': 49,
    'num_false_positives': 11,
    'num_true_negatives': 9,
    'num_false_negatives': 0,
    'precision': 0.8166666666666667,
    'recall': 1.0,
    'f1': 0.8990825688073394,
    'accuracy': None}