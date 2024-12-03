"""
Classifies: CHEBI:58946 acyl-CoA oxoanion
"""
from rdkit import Chem

def is_acyl_CoA_oxoanion(smiles: str):
    """
    Determines if a molecule is an acyl-CoA oxoanion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA oxoanion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of Coenzyme A (CoA) structure
    coa_substructure = Chem.MolFromSmarts("CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)([O-])=O)n1cnc2c(N)ncnc12)")
    if not mol.HasSubstructMatch(coa_substructure):
        return False, "Molecule does not contain Coenzyme A (CoA) structure"

    # Check for the presence of thioester bond (C(=O)S)
    thioester_substructure = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_substructure):
        return False, "Molecule does not contain thioester bond (C(=O)S)"

    # Check for deprotonated phosphate and/or diphosphate groups
    phosphate_substructure = Chem.MolFromSmarts("P(=O)([O-])([O-])")
    if not mol.HasSubstructMatch(phosphate_substructure):
        return False, "Molecule does not contain deprotonated phosphate groups"

    return True, "Molecule is an acyl-CoA oxoanion"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:58946',
                          'name': 'acyl-CoA oxoanion',
                          'definition': 'Any acyl coenzyme A thioester in '
                                        'which one or more of the phosphate '
                                        'and/or diphosphate groups has been '
                                        'deprotonated.',
                          'parents': ['CHEBI:58945']},
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
    'num_true_positives': 88,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}