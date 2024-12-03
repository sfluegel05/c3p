"""
Classifies: CHEBI:33247 organic group
"""
from rdkit import Chem

def is_organic_group(smiles: str):
    """
    Determines if a molecule is an organic group (any substituent group or skeleton containing carbon).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic group, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains at least one carbon atom
    if not any(atom.GetSymbol() == 'C' for atom in mol.GetAtoms()):
        return False, "No carbon atoms found in the molecule"

    return True, "Molecule contains carbon atoms, hence it is an organic group"

# Example usage
# print(is_organic_group("C([C@H](CO)N*)(=O)*"))  # L-serine residue
# print(is_organic_group("C(CCCCCCCCCCCCCCCCCCCC)CCCC*"))  # pentacosyl group


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33247',
                          'name': 'organic group',
                          'definition': 'Any substituent group or skeleton '
                                        'containing carbon.',
                          'parents': ['CHEBI:24433']},
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
    'num_true_positives': 165,
    'num_false_positives': 19,
    'num_true_negatives': 1,
    'num_false_negatives': 1,
    'precision': 0.8967391304347826,
    'recall': 0.9939759036144579,
    'f1': 0.9428571428571428,
    'accuracy': None}