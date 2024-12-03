"""
Classifies: CHEBI:60924 keratan sulfate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_keratan_sulfate(smiles: str):
    """
    Determines if a molecule is a keratan sulfate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a keratan sulfate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the repeating disaccharide unit of keratan sulfate
    keratan_sulfate_unit = Chem.MolFromSmarts("[C@H]1([C@H](O[C@H]([C@@H](CO)O[C@H]1O)CO)O)O")
    if keratan_sulfate_unit is None:
        return False, "Error in defining keratan sulfate unit"

    # Check for the presence of the repeating unit
    if not mol.HasSubstructMatch(keratan_sulfate_unit):
        return False, "No keratan sulfate repeating unit found"

    # Check for sulfo groups
    sulfo_group = Chem.MolFromSmarts("S(=O)(=O)[O-]")
    if sulfo_group is None:
        return False, "Error in defining sulfo group"

    if not mol.HasSubstructMatch(sulfo_group):
        return False, "No sulfo groups found"

    return True, "Keratan sulfate structure identified"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:60924',
                          'name': 'keratan sulfate',
                          'definition': 'A sulfated glycosaminoglycan, a '
                                        'linear polymer that consists of the '
                                        'repeating disaccharide  '
                                        '[3)-beta-Gal-(1->4)-beta-GlcNAc-(1->] '
                                        'and containing sulfo groups located '
                                        'at random positions.',
                          'parents': ['CHEBI:35722', 'CHEBI:37395']},
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 29,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}