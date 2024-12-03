"""
Classifies: CHEBI:24689 hydroxycinnamic acid
"""
from rdkit import Chem

def is_hydroxycinnamic_acid(smiles: str):
    """
    Determines if a molecule is a hydroxycinnamic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxycinnamic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the cinnamic acid core structure: C=CC(O)=O
    pattern = Chem.MolFromSmarts('C=CC(O)=O')
    if not mol.HasSubstructMatch(pattern):
        return False, "No cinnamic acid core structure found"

    # Check for hydroxy substituents on the aromatic ring
    hydroxy_pattern = Chem.MolFromSmarts('c[OH]')
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No hydroxy substituents found on the aromatic ring"

    return True, "Hydroxycinnamic acid identified"

# Example usage:
# result, reason = is_hydroxycinnamic_acid('CC(C)=CCc1cc(\C=C\C(O)=O)cc(C\C=C(/C)CO)c1O')
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24689',
                          'name': 'hydroxycinnamic acid',
                          'definition': 'Any member of the class of  cinnamic '
                                        'acids carrying one or more hydroxy '
                                        'substituents.',
                          'parents': ['CHEBI:23252']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[20:39:30] SMILES Parse Error: syntax error while parsing: '
             'O=C\x01C(CC(C/C1=C\\CC(CC(=O)N)CC(O)=O)C)C\n'
             '[20:39:30] SMILES Parse Error: Failed parsing SMILES '
             "'O=C\x01C(CC(C/C1=C\\CC(CC(=O)N)CC(O)=O)C)C' for input: "
             "'O=C\x01C(CC(C/C1=C\\CC(CC(=O)N)CC(O)=O)C)C'\n",
    'stdout': '',
    'num_true_positives': 18,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 10,
    'precision': 0.9,
    'recall': 0.6428571428571429,
    'f1': 0.75,
    'accuracy': None}