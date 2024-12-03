"""
Classifies: CHEBI:51737 alpha,beta-unsaturated carboxylic ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_alpha_beta_unsaturated_carboxylic_ester(smiles: str):
    """
    Determines if a molecule is an alpha,beta-unsaturated carboxylic ester.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha,beta-unsaturated carboxylic ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for alpha,beta-unsaturated carboxylic ester
    pattern1 = Chem.MolFromSmarts("C=CC(=O)OC")
    pattern2 = Chem.MolFromSmarts("C#CC(=O)OC")
    
    if mol.HasSubstructMatch(pattern1):
        return True, "Contains alpha,beta-unsaturated carboxylic ester pattern 1"
    elif mol.HasSubstructMatch(pattern2):
        return True, "Contains alpha,beta-unsaturated carboxylic ester pattern 2"
    else:
        return False, "Does not contain alpha,beta-unsaturated carboxylic ester pattern"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51737',
                          'name': 'alpha,beta-unsaturated carboxylic ester',
                          'definition': 'A carboxylic ester of general formula '
                                        'R(1)R(2)C=CR(3)-C(=O)OR(4) (R(4) =/= '
                                        'H) or R(1)C#C-C(=O)OR(2) (R(2) =/= H) '
                                        'in which the ester C=O function is '
                                        'conjugated to an unsaturated C-C bond '
                                        'at the alpha,beta position.',
                          'parents': ['CHEBI:33308']},
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
    'num_true_positives': 44,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9777777777777777,
    'f1': 0.9887640449438202,
    'accuracy': None}