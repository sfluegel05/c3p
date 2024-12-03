"""
Classifies: CHEBI:63551 carbohydrate acid derivative anion
"""
from rdkit import Chem

def is_carbohydrate_acid_derivative_anion(smiles: str):
    """
    Determines if a molecule is a carbohydrate acid derivative anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbohydrate acid derivative anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylate group (deprotonated carboxy group)
    carboxylate = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3:
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if neighbors.count('O') == 2 and neighbors.count('C') == 1:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O' and neighbor.GetFormalCharge() == -1:
                        carboxylate = True
                        break

    if not carboxylate:
        return False, "No carboxylate group found"

    # Check for the presence of multiple hydroxyl groups (common in carbohydrates)
    hydroxyl_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O':
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if neighbors.count('C') == 1:
                hydroxyl_count += 1

    if hydroxyl_count < 3:
        return False, "Not enough hydroxyl groups to be a carbohydrate derivative"

    # Check for glycosidic bonds (O-C-O pattern)
    glycosidic_bond = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O':
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if neighbors.count('C') == 2:
                glycosidic_bond = True
                break

    if not glycosidic_bond:
        return False, "No glycosidic bond found"

    return True, "Carbohydrate acid derivative anion"

# Example usage
smiles = "O([C@@H]1[C@H]([C@H](O[C@@H]2[C@H](O[C@@H](OC[C@@H]([C@@H](/C=C/CCCCCCCCCCCCC)O)NC(=O)*)[C@@H]([C@H]2O)O)CO)O[C@@H]([C@@H]1O)CO)O)[C@]3(O[C@]([C@@H]([C@H](C3)O)NC(C)=O)([C@@H]([C@H](O[C@]4(O[C@]([C@@H]([C@H](C4)O)NC(C)=O)([C@@H]([C@@H](CO)O)OC(C)=O)[H])C(=O)[O-])CO)O)[H])C([O-])=O"
result, reason = is_carbohydrate_acid_derivative_anion(smiles)
print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:63551',
                          'name': 'carbohydrate acid derivative anion',
                          'definition': 'A carboxylic acid anion resulting '
                                        'from the deprotonation of the carboxy '
                                        'group of a carbohydrate acid '
                                        'derivative.',
                          'parents': ['CHEBI:29067']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'True Carbohydrate acid derivative anion\n',
    'num_true_positives': 56,
    'num_false_positives': 3,
    'num_true_negatives': 17,
    'num_false_negatives': 4,
    'precision': 0.9491525423728814,
    'recall': 0.9333333333333333,
    'f1': 0.9411764705882353,
    'accuracy': None}