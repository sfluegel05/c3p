"""
Classifies: CHEBI:58944 dialkyl phosphate anion
"""
from rdkit import Chem

def is_dialkyl_phosphate_anion(smiles: str):
    """
    Determines if a molecule is a dialkyl phosphate anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dialkyl phosphate anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    phosphate_groups = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'P':
            neighbors = [n for n in atom.GetNeighbors()]
            if len(neighbors) == 4:
                o_count = sum(1 for n in neighbors if n.GetSymbol() == 'O')
                neg_o_count = sum(1 for n in neighbors if n.GetSymbol() == 'O' and n.GetFormalCharge() == -1)
                if o_count == 4 and neg_o_count >= 1:
                    phosphate_groups.append(atom)

    if not phosphate_groups:
        return False, "No phosphate group with appropriate oxygens found"

    for phosphate in phosphate_groups:
        neighbors = [n for n in phosphate.GetNeighbors() if n.GetSymbol() == 'O' and n.GetFormalCharge() == 0]
        if len(neighbors) < 2:
            continue

        alkyl_groups = []
        for oxygen in neighbors:
            for neighbor in oxygen.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    alkyl_groups.append(neighbor)

        if len(alkyl_groups) >= 2:
            return True, "Dialkyl phosphate anion found"

    return False, "No dialkyl phosphate anion structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:58944',
                          'name': 'dialkyl phosphate anion',
                          'definition': 'The conjugate base of a dialkyl '
                                        'phosphate compound',
                          'parents': ['CHEBI:58945']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 13,
    'num_false_positives': 3,
    'num_true_negatives': 10,
    'num_false_negatives': 0,
    'precision': 0.8125,
    'recall': 1.0,
    'f1': 0.896551724137931,
    'accuracy': None}