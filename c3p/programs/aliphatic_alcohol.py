"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol (an alcohol derived from an aliphatic compound).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of an alcohol group (hydroxyl group attached to a carbon)
    alcohol_groups = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and any(neigh.GetSymbol() == 'C' for neigh in atom.GetNeighbors())]
    if not alcohol_groups:
        return False, "No alcohol groups found"

    # Check if the molecule is aliphatic (contains no aromatic rings)
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Molecule contains aromatic rings"

    # Ensure the alcohol group is attached to an aliphatic carbon
    for alcohol in alcohol_groups:
        for neighbor in alcohol.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and not neighbor.GetIsAromatic():
                return True, "Molecule is an aliphatic alcohol"

    return False, "Alcohol group is not attached to an aliphatic carbon"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:2571',
                          'name': 'aliphatic alcohol',
                          'definition': 'An  alcohol derived from an aliphatic '
                                        'compound.',
                          'parents': ['CHEBI:30879']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '[20:51:44] SMILES Parse Error: unclosed ring for input: '
             "'[H][C@]12OC(=O)C(=C)[C@]1([H])[C@H](O)C\\C(C)=C\\CC\\C(C)=C\x02'\n",
    'stdout': '',
    'num_true_positives': 80,
    'num_false_positives': 13,
    'num_true_negatives': 7,
    'num_false_negatives': 10,
    'precision': 0.8602150537634409,
    'recall': 0.8888888888888888,
    'f1': 0.8743169398907102,
    'accuracy': None}