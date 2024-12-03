"""
Classifies: CHEBI:60971 aminophospholipid
"""
from rdkit import Chem

def is_aminophospholipid(smiles: str):
    """
    Determines if a molecule is an aminophospholipid (a phospholipid that contains one or more amino groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aminophospholipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a phosphate group
    phosphate_group = Chem.MolFromSmarts('P(=O)(O)(O)O')
    if not mol.HasSubstructMatch(phosphate_group):
        return False, "No phosphate group found"

    # Check for the presence of at least one amino group
    amino_group = Chem.MolFromSmarts('[NX3;H2,H1;!$(NC=O)]')
    if not mol.HasSubstructMatch(amino_group):
        return False, "No amino group found"

    return True, "Molecule is an aminophospholipid"

# Example usage:
# smiles = "P(OCC(OC(=O)CCCCCCC/C=C\CCCCCCCC)COC(=O)CCCCCCCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O"
# print(is_aminophospholipid(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:60971',
                          'name': 'aminophospholipid',
                          'definition': 'A phospholipid that contains one or '
                                        'more amino groups.',
                          'parents': ['CHEBI:16247', 'CHEBI:50047']},
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
    'num_true_positives': 44,
    'num_false_positives': 6,
    'num_true_negatives': 14,
    'num_false_negatives': 0,
    'precision': 0.88,
    'recall': 1.0,
    'f1': 0.9361702127659575,
    'accuracy': None}