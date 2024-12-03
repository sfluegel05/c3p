"""
Classifies: CHEBI:67265 beta-diketone
"""
from rdkit import Chem

def is_beta_diketone(smiles: str):
    """
    Determines if a molecule is a beta-diketone (a diketone in which the two keto groups are separated by a single carbon atom).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-diketone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all carbonyl groups (C=O)
    carbonyls = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble() == 2:
                    carbonyls.append(atom.GetIdx())
                    break

    if len(carbonyls) < 2:
        return False, "Less than two carbonyl groups found"

    # Check if there are two carbonyl groups separated by a single carbon atom
    for i in range(len(carbonyls)):
        for j in range(i + 1, len(carbonyls)):
            path = Chem.rdmolops.GetShortestPath(mol, carbonyls[i], carbonyls[j])
            if len(path) == 3:
                # Check if the central atom is a carbon
                if mol.GetAtomWithIdx(path[1]).GetSymbol() == 'C':
                    return True, "Beta-diketone structure found"

    return False, "No beta-diketone structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:67265',
                          'name': 'beta-diketone',
                          'definition': 'A diketone in which the two keto '
                                        'groups are separated by a single '
                                        'carbon atom.',
                          'parents': ['CHEBI:46640']},
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
    'num_true_positives': 13,
    'num_false_positives': 1,
    'num_true_negatives': 12,
    'num_false_negatives': 0,
    'precision': 0.9285714285714286,
    'recall': 1.0,
    'f1': 0.962962962962963,
    'accuracy': None}