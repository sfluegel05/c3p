"""
Classifies: CHEBI:36709 aminoquinoline
"""
from rdkit import Chem

def is_aminoquinoline(smiles: str):
    """
    Determines if a molecule is an aminoquinoline (quinoline skeleton substituted by one or more amino or substituted-amino groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aminoquinoline, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for quinoline skeleton
    quinoline = Chem.MolFromSmarts('c1ccc2ncccc2c1')
    if not mol.HasSubstructMatch(quinoline):
        return False, "No quinoline skeleton found"

    # Check for amino groups attached to the quinoline skeleton
    amino_group = Chem.MolFromSmarts('[NX3][CX4]')
    substituted_amino_group = Chem.MolFromSmarts('[NX3;!$(NC=O)]')

    quinoline_matches = mol.GetSubstructMatches(quinoline)
    for match in quinoline_matches:
        quinoline_atoms = set(match)

        for atom_idx in quinoline_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in quinoline_atoms:
                    if neighbor.GetSmarts() == 'N' or neighbor.GetSmarts() == '[NX3;!$(NC=O)]':
                        return True, "Aminoquinoline found"

    return False, "No amino or substituted-amino groups found attached to the quinoline skeleton"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36709',
                          'name': 'aminoquinoline',
                          'definition': 'Any member of the class of quinolines '
                                        'in which the quinoline skeleton is '
                                        'substituted by one or more amino or '
                                        'substituted-amino groups.',
                          'parents': ['CHEBI:26513', 'CHEBI:33860']},
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
    'num_true_positives': 4,
    'num_false_positives': 0,
    'num_true_negatives': 17,
    'num_false_negatives': 13,
    'precision': 1.0,
    'recall': 0.23529411764705882,
    'f1': 0.38095238095238093,
    'accuracy': None}