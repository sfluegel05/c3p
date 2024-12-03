"""
Classifies: CHEBI:78189 alkylacylglycero-3-phosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_alkylacylglycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is an alkylacylglycero-3-phosphocholine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkylacylglycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the glycerophosphocholine core structure
    core_smarts = '[C@H](COP([O-])(=O)OCC[N+](C)(C)C)O'
    core = Chem.MolFromSmarts(core_smarts)
    if not mol.HasSubstructMatch(core):
        return False, "Does not contain glycerophosphocholine core structure"

    # Check for the presence of one alkyl group (simple alkane chain)
    alkyl_smarts = 'OCC[*]'
    alkyl = Chem.MolFromSmarts(alkyl_smarts)
    if not mol.HasSubstructMatch(alkyl):
        return False, "Does not contain alkyl group"

    # Check for the presence of one acyl group (ester linkage to fatty acid chain)
    acyl_smarts = 'OC(=O)[*]'
    acyl = Chem.MolFromSmarts(acyl_smarts)
    if not mol.HasSubstructMatch(acyl):
        return False, "Does not contain acyl group"

    return True, "Contains glycerophosphocholine core with one alkyl and one acyl group in unspecified positions"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:78189',
                          'name': 'alkylacylglycero-3-phosphocholine',
                          'definition': 'A glycerophosphocholine that is '
                                        'glycero-3-phosphocholine carrying one '
                                        'alkyl and one acyl group in '
                                        'unspecified positions. Major species '
                                        'at pH 7.3.',
                          'parents': ['CHEBI:36313', 'CHEBI:64611']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 18,
    'num_false_positives': 20,
    'num_true_negatives': 0,
    'num_false_negatives': 8,
    'precision': 0.47368421052631576,
    'recall': 0.6923076923076923,
    'f1': 0.5625,
    'accuracy': None}