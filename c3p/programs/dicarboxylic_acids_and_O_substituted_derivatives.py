"""
Classifies: CHEBI:131927 dicarboxylic acids and O-substituted derivatives
"""
from rdkit import Chem

def is_dicarboxylic_acids_and_O_substituted_derivatives(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid or an O-substituted derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule fits the class, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    carboxyl_groups = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if neighbors.count('O') == 2 and neighbors.count('C') == 1:
                carboxyl_groups.append(atom)

    if len(carboxyl_groups) < 2:
        return False, "Less than two carboxyl groups found"

    for atom in carboxyl_groups:
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'O':
                for sub_neighbor in neighbor.GetNeighbors():
                    if sub_neighbor.GetSymbol() != 'C' and sub_neighbor.GetSymbol() != 'H':
                        return True, f"O-substituted derivative with substituent: {sub_neighbor.GetSymbol()}"

    return True, "Dicarboxylic acid without O-substitution"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:131927',
                          'name': 'dicarboxylic acids and O-substituted '
                                  'derivatives',
                          'definition': 'A class of carbonyl compound '
                                        'encompassing dicarboxylic acids and '
                                        'any derivatives obtained by '
                                        'substitution of either one or both of '
                                        'the carboxy hydrogens.',
                          'parents': ['CHEBI:36586']},
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
    'num_true_positives': 91,
    'num_false_positives': 4,
    'num_true_negatives': 16,
    'num_false_negatives': 0,
    'precision': 0.9578947368421052,
    'recall': 1.0,
    'f1': 0.978494623655914,
    'accuracy': None}