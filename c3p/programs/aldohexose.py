"""
Classifies: CHEBI:33917 aldohexose
"""
from rdkit import Chem

def is_aldohexose(smiles: str):
    """
    Determines if a molecule is an aldohexose (a hexose with a potential aldehyde group at one end).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldohexose, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 6 carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if num_carbons != 6:
        return False, "Molecule does not have 6 carbon atoms"

    # Check for an aldehyde group
    aldehyde_group = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3:
            neighbors = atom.GetNeighbors()
            if any(neighbor.GetSymbol() == 'O' and neighbor.GetTotalDegree() == 1 for neighbor in neighbors):
                aldehyde_group = True
                break

    if not aldehyde_group:
        return False, "No aldehyde group found"

    # Check for hydroxyl groups
    hydroxyl_groups = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and any(neighbor.GetSymbol() == 'C' for neighbor in atom.GetNeighbors()))
    if hydroxyl_groups < 4:
        return False, "Molecule does not have sufficient hydroxyl groups"

    return True, "Molecule is an aldohexose"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33917',
                          'name': 'aldohexose',
                          'definition': 'A hexose with a (potential) aldehyde '
                                        'group at one end.',
                          'parents': ['CHEBI:15693', 'CHEBI:18133']},
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
    'num_false_positives': 1,
    'num_true_negatives': 15,
    'num_false_negatives': 12,
    'precision': 0.8,
    'recall': 0.25,
    'f1': 0.38095238095238093,
    'accuracy': None}