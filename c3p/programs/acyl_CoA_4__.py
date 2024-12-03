"""
Classifies: CHEBI:58342 acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is an acyl-CoA(4-) (An acyl-CoA oxoanion arising from deprotonation of the phosphate and diphosphate OH groups of any acyl-CoA; major species at pH 7.3).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the CoA core structure
    pattern = Chem.MolFromSmarts("C(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)([O-])[O-])n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)")
    if not mol.HasSubstructMatch(pattern):
        return False, "Missing CoA core structure"

    # Check for the presence of the acyl group
    acyl_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)[C@H](O)")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "Missing acyl group"

    # Check for the presence of the oxoanion groups (deprotonated phosphate and diphosphate)
    oxoanion_pattern = Chem.MolFromSmarts("P(=O)([O-])OP(=O)([O-])")
    if not mol.HasSubstructMatch(oxoanion_pattern):
        return False, "Missing oxoanion groups"

    return True, "Valid acyl-CoA(4-) structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:58342',
                          'name': 'acyl-CoA(4-)',
                          'definition': 'An acyl-CoA oxoanion arising from '
                                        'deprotonation of the phosphate and '
                                        'diphosphate OH groups of any '
                                        'acyl-CoA; major species at pH 7.3.',
                          'parents': ['CHEBI:58946']},
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
    'num_true_positives': 77,
    'num_false_positives': 3,
    'num_true_negatives': 17,
    'num_false_negatives': 1,
    'precision': 0.9625,
    'recall': 0.9871794871794872,
    'f1': 0.9746835443037976,
    'accuracy': None}