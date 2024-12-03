"""
Classifies: CHEBI:61902 hydroxy fatty acyl-CoA
"""
from rdkit import Chem

def is_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acyl-CoA.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the CoA moiety
    coa_pattern = Chem.MolFromSmarts('C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12')
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Check for the presence of a hydroxy fatty acid moiety
    hydroxy_fatty_acid_pattern = Chem.MolFromSmarts('C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C')
    if not mol.HasSubstructMatch(hydroxy_fatty_acid_pattern):
        return False, "Hydroxy fatty acid moiety not found"

    # Check for the presence of a hydroxy group in the fatty acid chain
    hydroxy_group_pattern = Chem.MolFromSmarts('[C@H](O)')
    if not mol.HasSubstructMatch(hydroxy_group_pattern):
        return False, "Hydroxy group not found in fatty acid chain"

    return True, "Molecule is a hydroxy fatty acyl-CoA"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61902',
                          'name': 'hydroxy fatty acyl-CoA',
                          'definition': 'A fatty-acyl-CoA that results from '
                                        'the formal condensation of the thiol '
                                        'group of coenzyme A with the carboxy '
                                        'group of any hydroxy fatty acid.',
                          'parents': ['CHEBI:37554', 'CHEBI:62618']},
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
    'num_true_positives': 12,
    'num_false_positives': 12,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.5,
    'recall': 1.0,
    'f1': 0.6666666666666666,
    'accuracy': None}