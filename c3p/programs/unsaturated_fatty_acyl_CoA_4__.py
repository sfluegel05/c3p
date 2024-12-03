"""
Classifies: CHEBI:77632 unsaturated fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_unsaturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA(4-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the CoA structure
    coa_smarts = "C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    coa = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa):
        return False, "No CoA structure found"

    # Check for the presence of a fatty acyl chain
    fatty_acyl_smarts = "C(=O)SCCNC(=O)CCNC(=O)[C@H](O)"
    fatty_acyl = Chem.MolFromSmarts(fatty_acyl_smarts)
    if not mol.HasSubstructMatch(fatty_acyl):
        return False, "No fatty acyl chain found"

    # Check for unsaturation (presence of double bonds in the fatty acyl chain)
    unsaturation_smarts = "C=C"
    unsaturation = Chem.MolFromSmarts(unsaturation_smarts)
    if not mol.HasSubstructMatch(unsaturation):
        return False, "No double bonds found in the fatty acyl chain"

    # Check for deprotonation of the phosphate and diphosphate OH groups
    deprotonation_smarts = "P([O-])(=O)OP([O-])(=O)"
    deprotonation = Chem.MolFromSmarts(deprotonation_smarts)
    if not mol.HasSubstructMatch(deprotonation):
        return False, "No deprotonation of phosphate and diphosphate OH groups found"

    return True, "Molecule is an unsaturated fatty acyl-CoA(4-)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:77632',
                          'name': 'unsaturated fatty acyl-CoA(4-)',
                          'definition': 'A fatty acyl-CoA(4-) arising from '
                                        'deprotonation of the phosphate and '
                                        'diphosphate OH groups of any '
                                        'unsaturated fatty acyl-CoA; major '
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
    'num_true_positives': 14,
    'num_false_positives': 2,
    'num_true_negatives': 15,
    'num_false_negatives': 3,
    'precision': 0.875,
    'recall': 0.8235294117647058,
    'f1': 0.8484848484848485,
    'accuracy': None}