"""
Classifies: CHEBI:77636 fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA(4-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the acyl-CoA core structure
    acyl_CoA_pattern = Chem.MolFromSmarts('C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12')
    if not mol.HasSubstructMatch(acyl_CoA_pattern):
        return False, "Does not contain the acyl-CoA core structure"

    # Check for the presence of fatty acyl chain
    fatty_acyl_chain_pattern = Chem.MolFromSmarts('C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12')
    if not mol.HasSubstructMatch(fatty_acyl_chain_pattern):
        return False, "Does not contain a fatty acyl chain"

    # Check for deprotonation of phosphate and diphosphate OH groups
    phosphate_pattern = Chem.MolFromSmarts('OP([O-])(=O)')
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group is not deprotonated"
    
    diphosphate_pattern = Chem.MolFromSmarts('OP([O-])(=O)OP([O-])(=O)')
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "Diphosphate group is not deprotonated"

    return True, "Molecule is a fatty acyl-CoA(4-)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:77636',
                          'name': 'fatty acyl-CoA(4-)',
                          'definition': 'An acyl-CoA(4-) arising from '
                                        'deprotonation of the phosphate and '
                                        'diphosphate OH groups of any fatty '
                                        'acyl-CoA; major species at pH 7.3.',
                          'parents': ['CHEBI:58342']},
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
    'num_true_positives': 47,
    'num_false_positives': 11,
    'num_true_negatives': 9,
    'num_false_negatives': 0,
    'precision': 0.8103448275862069,
    'recall': 1.0,
    'f1': 0.8952380952380952,
    'accuracy': None}