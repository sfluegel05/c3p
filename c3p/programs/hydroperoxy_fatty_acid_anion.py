"""
Classifies: CHEBI:64012 hydroperoxy fatty acid anion
"""
from rdkit import Chem

def is_hydroperoxy_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a hydroperoxy fatty acid anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroperoxy fatty acid anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the carboxylate anion group [C(=O)[O-]]
    carboxylate_pattern = Chem.MolFromSmarts('C(=O)[O-]')
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate anion group found"

    # Check for the hydroperoxy group [COO]
    hydroperoxy_pattern = Chem.MolFromSmarts('COO')
    if not mol.HasSubstructMatch(hydroperoxy_pattern):
        return False, "No hydroperoxy group found"

    # Check for the presence of double bonds (polyunsaturated)
    double_bond_pattern = Chem.MolFromSmarts('C=C')
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No double bonds found"

    # Ensure the molecule is a fatty acid (long carbon chain)
    carbon_chain_length = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            carbon_chain_length += 1
    if carbon_chain_length < 8:  # typically fatty acids have at least 8 carbons
        return False, "Carbon chain is too short to be a fatty acid"

    return True, "Molecule is a hydroperoxy fatty acid anion"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:64012',
                          'name': 'hydroperoxy fatty acid anion',
                          'definition': 'A fatty acid anion that is the '
                                        'conjugate base of any hydroperoxy '
                                        'fatty acid, formed by deprotonation '
                                        'of the carboxylic acid moiety.',
                          'parents': ['CHEBI:28868']},
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
    'num_true_positives': 9,
    'num_false_positives': 0,
    'num_true_negatives': 10,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9,
    'f1': 0.9473684210526316,
    'accuracy': None}