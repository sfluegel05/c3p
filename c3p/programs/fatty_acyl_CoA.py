"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for Coenzyme A moiety
    coA_smiles = "CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCS"
    coA_mol = Chem.MolFromSmiles(coA_smiles)
    if coA_mol is None:
        return False, "Error in Coenzyme A SMILES string"

    # Check if the molecule contains the Coenzyme A moiety
    if not mol.HasSubstructMatch(coA_mol):
        return False, "Coenzyme A moiety not found"

    # Check for fatty acid moiety
    fatty_acid_smiles = "C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    fatty_acid_mol = Chem.MolFromSmiles(fatty_acid_smiles)
    if fatty_acid_mol is None:
        return False, "Error in fatty acid SMILES string"

    # Check if the molecule contains the fatty acid moiety
    if not mol.HasSubstructMatch(fatty_acid_mol):
        return False, "Fatty acid moiety not found"

    return True, "Molecule is a fatty acyl-CoA"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37554',
                          'name': 'fatty acyl-CoA',
                          'definition': 'An acyl-CoA that results from the '
                                        'formal condensation of the thiol '
                                        'group of coenzyme A with the carboxy '
                                        'group of any fatty acid.',
                          'parents': ['CHEBI:17984', 'CHEBI:61697']},
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
    'num_true_positives': 53,
    'num_false_positives': 8,
    'num_true_negatives': 12,
    'num_false_negatives': 0,
    'precision': 0.8688524590163934,
    'recall': 1.0,
    'f1': 0.9298245614035088,
    'accuracy': None}