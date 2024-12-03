"""
Classifies: CHEBI:138366 bile acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_bile_acids(smiles: str):
    """
    Determines if a molecule is a bile acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is a steroid
    steroid_core = Chem.MolFromSmarts('C1CCC2C3CCC4CC(C)(C)CCC4(C)C3CCC21')
    if not mol.HasSubstructMatch(steroid_core):
        return False, "Molecule is not a steroid"

    # Check for hydroxy groups
    hydroxy_groups = mol.GetSubstructMatches(Chem.MolFromSmarts('O'))
    if len(hydroxy_groups) < 1:
        return False, "Molecule does not contain hydroxy groups"

    # Check for carboxylic acid or its derivatives
    carboxylic_acid = Chem.MolFromSmarts('C(=O)O')
    carboxylic_acid_derivatives = Chem.MolFromSmarts('C(=O)[N,O]')
    if not (mol.HasSubstructMatch(carboxylic_acid) or mol.HasSubstructMatch(carboxylic_acid_derivatives)):
        return False, "Molecule does not contain carboxylic acid or its derivatives"

    # Check for 5beta or 5alpha configuration
    five_beta_config = Chem.MolFromSmarts('[C@@H]1CC[C@](C)(CC1)')
    five_alpha_config = Chem.MolFromSmarts('[C@H]1CC[C@@](C)(CC1)')
    if not (mol.HasSubstructMatch(five_beta_config) or mol.HasSubstructMatch(five_alpha_config)):
        return False, "Molecule does not have 5beta or 5alpha configuration"

    return True, "Molecule is a bile acid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:138366',
                          'name': 'bile acids',
                          'definition': 'Any member of a group of hydroxy '
                                        'steroids occuring in bile, where they '
                                        'are present as the sodium salts of '
                                        'their amides with glycine or taurine. '
                                        'In mammals bile acids almost '
                                        'invariably have 5beta-configuration, '
                                        'while in lower vertebrates, some bile '
                                        'acids, known as allo-bile acids, have '
                                        '5alpha-configuration.',
                          'parents': [   'CHEBI:25384',
                                         'CHEBI:35350',
                                         'CHEBI:36078']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[19:20:29] SMILES Parse Error: syntax error while parsing: '
             'C(CC/C=C\\C=C\x01/O[C@H]1C/C=C\\CC)CCCCC(=O)O\n'
             '[19:20:29] SMILES Parse Error: Failed parsing SMILES '
             "'C(CC/C=C\\C=C\x01/O[C@H]1C/C=C\\CC)CCCCC(=O)O' for input: "
             "'C(CC/C=C\\C=C\x01/O[C@H]1C/C=C\\CC)CCCCC(=O)O'\n",
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 42,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}