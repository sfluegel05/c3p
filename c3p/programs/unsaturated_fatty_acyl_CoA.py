"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for Coenzyme A (CoA) part
    coa_fragment = Chem.MolFromSmarts('C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12')
    if not mol.HasSubstructMatch(coa_fragment):
        return False, "No Coenzyme A (CoA) part found"

    # Check for unsaturated fatty acid part
    unsaturated_fatty_acid_fragment = Chem.MolFromSmarts('C(=O)SCCNC(=O)CCNC(=O)[C@H](O)')
    if not mol.HasSubstructMatch(unsaturated_fatty_acid_fragment):
        return False, "No unsaturated fatty acid part found"

    # Check for double bonds in the fatty acid chain
    double_bond_fragment = Chem.MolFromSmarts('C=C')
    if not mol.HasSubstructMatch(double_bond_fragment):
        return False, "No double bonds found in the fatty acid chain"

    return True, "Molecule is an unsaturated fatty acyl-CoA"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51006',
                          'name': 'unsaturated fatty acyl-CoA',
                          'definition': 'A fatty acyl-CoA that results from '
                                        'the formal condensation of the thiol '
                                        'group of coenzyme A with the carboxy '
                                        'group of any unsaturated fatty acid.',
                          'parents': ['CHEBI:37554']},
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
    'num_true_positives': 33,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 0,
    'precision': 0.9428571428571428,
    'recall': 1.0,
    'f1': 0.9705882352941176,
    'accuracy': None}