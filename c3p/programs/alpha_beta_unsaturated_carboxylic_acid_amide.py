"""
Classifies: CHEBI:51750 alpha,beta-unsaturated carboxylic acid amide
"""
from rdkit import Chem

def is_alpha_beta_unsaturated_carboxylic_acid_amide(smiles: str):
    """
    Determines if a molecule is an alpha,beta-unsaturated carboxylic acid amide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha,beta-unsaturated carboxylic acid amide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the amide group
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found"

    # Check for the conjugated C=C bond in alpha,beta position or C#C bond
    alpha_beta_unsat_pattern = Chem.MolFromSmarts("C=CC(=O)N")
    conjugated_triple_bond_pattern = Chem.MolFromSmarts("C#CC(=O)N")

    if mol.HasSubstructMatch(alpha_beta_unsat_pattern):
        return True, "Contains alpha,beta-unsaturated C=C bond conjugated to amide group"
    elif mol.HasSubstructMatch(conjugated_triple_bond_pattern):
        return True, "Contains conjugated C#C bond to amide group"

    return False, "No alpha,beta-unsaturated C=C or C#C bond conjugated to amide group found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51750',
                          'name': 'alpha,beta-unsaturated carboxylic acid '
                                  'amide',
                          'definition': 'A monocarboxylic amide of general '
                                        'formula '
                                        'R(1)R(2)C=CR(3)-C(=O)NR(4)R(5)  or '
                                        'R(1)C#C-C(=O)NR(2)R(3) in which the '
                                        'amide C=O function is conjugated to '
                                        'an unsaturated C-C bond at the '
                                        'alpha,beta position.',
                          'parents': ['CHEBI:37622']},
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
    'num_true_positives': 18,
    'num_false_positives': 0,
    'num_true_negatives': 18,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}