"""
Classifies: CHEBI:84465 acyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_acyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is an acyl-sn-glycero-3-phosphocholine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the glycerophosphocholine backbone
    pattern_gpc = Chem.MolFromSmarts('OCC[N+](C)(C)C')
    if not mol.HasSubstructMatch(pattern_gpc):
        return False, "No glycerophosphocholine backbone found"

    # Check for the presence of the phosphate group
    pattern_phosphate = Chem.MolFromSmarts('P(=O)(O)(O)')
    if not mol.HasSubstructMatch(pattern_phosphate):
        return False, "No phosphate group found"

    # Check for the presence of the glycerol backbone with stereochemistry
    pattern_glycerol = Chem.MolFromSmarts('O[C@H](CO)CO')
    if not mol.HasSubstructMatch(pattern_glycerol):
        return False, "No glycerol backbone with defined stereochemistry found"

    # Check for the acyl group attached to the glycerol backbone
    pattern_acyl = Chem.MolFromSmarts('O=C')
    if not mol.HasSubstructMatch(pattern_acyl):
        return False, "No acyl group found"

    return True, "Molecule is an acyl-sn-glycero-3-phosphocholine"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:84465',
                          'name': 'acyl-sn-glycero-3-phosphocholine',
                          'definition': 'A lysophosphatidylcholine with '
                                        'defined stereochemistry where the '
                                        'position of the acyl group is '
                                        'unknown.',
                          'parents': ['CHEBI:60479']},
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
    'num_false_positives': 4,
    'num_true_negatives': 7,
    'num_false_negatives': 0,
    'precision': 0.7333333333333333,
    'recall': 1.0,
    'f1': 0.846153846153846,
    'accuracy': None}