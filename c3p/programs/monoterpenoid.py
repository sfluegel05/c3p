"""
Classifies: CHEBI:25409 monoterpenoid
"""
from rdkit import Chem

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains 10 carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if not (8 <= num_carbons <= 10):
        return False, f"Number of carbon atoms is {num_carbons}, not in the range of 8-10"

    # Check if the molecule is derived from a monoterpene
    # This is a simplified check to see if the molecule has a C10 skeleton
    # and contains typical monoterpenoid functional groups (e.g., alcohols, ketones, aldehydes, carboxylic acids)
    functional_groups = ['O', 'N', 'S']  # Common functional groups in monoterpenoids
    contains_functional_groups = any(atom.GetSymbol() in functional_groups for atom in mol.GetAtoms())
    
    if contains_functional_groups:
        return True, "Molecule is a monoterpenoid with functional groups"

    return False, "Molecule does not contain typical monoterpenoid functional groups"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25409',
                          'name': 'monoterpenoid',
                          'definition': 'Any terpenoid derived from a '
                                        'monoterpene. The term includes '
                                        'compounds in which the C10 skeleton '
                                        'of the parent monoterpene has been '
                                        'rearranged or modified by the removal '
                                        'of one or more skeletal atoms '
                                        '(generally methyl groups).',
                          'parents': ['CHEBI:26873']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '[13:37:17] SMILES Parse Error: syntax error while parsing: '
             'Br/C=C\x01/C2=C3C(=C[C@]4(C[C@@H]([C@@]([C@H]4[C@@H]3CC1)(C[C@H](O)[C@H](O)[C@@](O)(CO)C)C)O)C)[C@H](C2)C\n'
             '[13:37:17] SMILES Parse Error: Failed parsing SMILES '
             "'Br/C=C\x01/C2=C3C(=C[C@]4(C[C@@H]([C@@]([C@H]4[C@@H]3CC1)(C[C@H](O)[C@H](O)[C@@](O)(CO)C)C)O)C)[C@H](C2)C' "
             'for input: '
             "'Br/C=C\x01/C2=C3C(=C[C@]4(C[C@@H]([C@@]([C@H]4[C@@H]3CC1)(C[C@H](O)[C@H](O)[C@@](O)(CO)C)C)O)C)[C@H](C2)C'\n",
    'stdout': '',
    'num_true_positives': 37,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 81,
    'precision': 1.0,
    'recall': 0.3135593220338983,
    'f1': 0.4774193548387097,
    'accuracy': None}