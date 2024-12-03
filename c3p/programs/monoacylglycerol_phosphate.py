"""
Classifies: CHEBI:16961 monoacylglycerol phosphate
"""
from rdkit import Chem

def is_monoacylglycerol_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol phosphate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacylglycerol phosphate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate group
    phosphate_group = Chem.MolFromSmarts('P(O)(O)=O')
    if not mol.HasSubstructMatch(phosphate_group):
        return False, "No phosphate group found"

    # Check for glycerol backbone
    glycerol_backbone = Chem.MolFromSmarts('OCC(O)CO')
    if not mol.HasSubstructMatch(glycerol_backbone):
        return False, "No glycerol backbone found"

    # Check for one fatty acid ester linkage
    ester_linkage = Chem.MolFromSmarts('C(=O)O')
    ester_count = len(mol.GetSubstructMatches(ester_linkage))
    if ester_count != 1:
        return False, f"Expected 1 ester linkage, found {ester_count}"

    return True, "Molecule is a monoacylglycerol phosphate"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16961',
                          'name': 'monoacylglycerol phosphate',
                          'definition': 'Derivatives of phosphoglycerols which '
                                        'have only one of the alcohol groups '
                                        'of the glycerol backbone ester-linked '
                                        'with a fatty acid.',
                          'parents': ['CHEBI:26707', 'CHEBI:37739']},
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
    'num_true_positives': 11,
    'num_false_positives': 5,
    'num_true_negatives': 8,
    'num_false_negatives': 2,
    'precision': 0.6875,
    'recall': 0.8461538461538461,
    'f1': 0.7586206896551724,
    'accuracy': None}