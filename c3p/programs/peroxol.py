"""
Classifies: CHEBI:35924 peroxol
"""
from rdkit import Chem

def is_peroxol(smiles: str):
    """
    Determines if a molecule is a peroxol (monosubstitution product of hydrogen peroxide HOOH, having the skeleton ROOH, in which R is any organyl group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a peroxol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the hydroperoxy group (ROOH)
    found_hydroperoxy = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O':
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 2:
                if neighbors[0].GetSymbol() == 'O' and neighbors[1].GetSymbol() == 'H':
                    found_hydroperoxy = True
                    break
                elif neighbors[1].GetSymbol() == 'O' and neighbors[0].GetSymbol() == 'H':
                    found_hydroperoxy = True
                    break

    if not found_hydroperoxy:
        return False, "No hydroperoxy group (ROOH) found"

    # Check if the other part is an organyl group (organic group)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O':
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 2:
                if neighbors[0].GetSymbol() == 'O' and neighbors[1].GetSymbol() == 'H':
                    organyl = neighbors[0].GetNeighbors()
                    if any(neighbor.GetSymbol() != 'H' for neighbor in organyl):
                        return True, "Valid peroxol with organyl group"
                elif neighbors[1].GetSymbol() == 'O' and neighbors[0].GetSymbol() == 'H':
                    organyl = neighbors[1].GetNeighbors()
                    if any(neighbor.GetSymbol() != 'H' for neighbor in organyl):
                        return True, "Valid peroxol with organyl group"

    return False, "No valid organyl group found"

# Example usage:
smiles = "CC(C)(C)OO"  # tert-butyl hydroperoxide
print(is_peroxol(smiles))  # Should return (True, "Valid peroxol with organyl group")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35924',
                          'name': 'peroxol',
                          'definition': 'Monosubstitution products of hydrogen '
                                        'peroxide HOOH, having the skeleton '
                                        'ROOH, in which R is any organyl '
                                        'group.',
                          'parents': [   'CHEBI:35923',
                                         'CHEBI:36963',
                                         'CHEBI:37863']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[23:43:27] Explicit valence for atom # 1 Si, 8, is greater than '
             'permitted\n',
    'stdout': "(False, 'No hydroperoxy group (ROOH) found')\n",
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 17,
    'num_false_negatives': 17,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}