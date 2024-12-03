"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the CoA moiety
    coa_smarts = "C(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)([O-])[O-])n2cnc3c(N)ncnc23)"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Check for the presence of the fatty acyl chain
    fatty_acyl_smarts = "C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)"
    fatty_acyl_pattern = Chem.MolFromSmarts(fatty_acyl_smarts)
    if not mol.HasSubstructMatch(fatty_acyl_pattern):
        return False, "No fatty acyl chain found"

    # Check for the deprotonation state (4- charge)
    phosphate_smarts = "P(=O)([O-])[O-]"
    phosphate_pattern = Chem.MolFromSmarts(phosphate_smarts)
    if len(mol.GetSubstructMatches(phosphate_pattern)) < 2:
        return False, "Phosphate groups not deprotonated"

    diphosphate_smarts = "P(=O)([O-])OP(=O)([O-])[O-]"
    diphosphate_pattern = Chem.MolFromSmarts(diphosphate_smarts)
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "Diphosphate groups not deprotonated"

    return True, "Valid long-chain fatty acyl-CoA(4-)"

# Example usage:
# result, reason = is_long_chain_fatty_acyl_CoA_4__("CCCCC\C=C/C\C=C/CCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12")
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:83139',
                          'name': 'long-chain fatty acyl-CoA(4-)',
                          'definition': 'A fatty acyl-CoA(4-) arising from '
                                        'deprotonation of the phosphate and '
                                        'diphosphate OH groups of any '
                                        'long-chain fatty acyl-CoA; major '
                                        'species at pH 7.3.',
                          'parents': ['CHEBI:77636']},
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
    'num_true_negatives': 17,
    'num_false_negatives': 17,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}