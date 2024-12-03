"""
Classifies: CHEBI:62618 hydroxyacyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_hydroxyacyl_CoA(smiles: str):
    """
    Determines if a molecule is a hydroxyacyl-CoA.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxyacyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains a Coenzyme A moiety
    coA_smarts = Chem.MolFromSmarts('CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC')
    if not mol.HasSubstructMatch(coA_smarts):
        return False, "No Coenzyme A moiety found"

    # Check if the molecule contains a hydroxy group attached to an acyl group
    hydroxyacyl_smarts = Chem.MolFromSmarts('[CX3](=O)[CX4][OX2H]')
    if not mol.HasSubstructMatch(hydroxyacyl_smarts):
        return False, "No hydroxyacyl group found"

    return True, "Molecule is a hydroxyacyl-CoA"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:62618',
                          'name': 'hydroxyacyl-CoA',
                          'definition': 'An acyl-CoA that results from the '
                                        'formal condensation of the thiol '
                                        'group of coenzyme A with the carboxy '
                                        'group of any hydroxycarboxylic acid.',
                          'parents': ['CHEBI:17984']},
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
    'num_true_positives': 17,
    'num_false_positives': 17,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.5,
    'recall': 1.0,
    'f1': 0.6666666666666666,
    'accuracy': None}