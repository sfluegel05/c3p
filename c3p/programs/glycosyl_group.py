"""
Classifies: CHEBI:24403 glycosyl group
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glycosyl_group(smiles: str):
    """
    Determines if a molecule is a glycosyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosyl group, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a hemiacetal carbon (carbon connected to two oxygens, one of which is an ether and the other an alcohol)
    hemiacetal_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            oxygen_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'O']
            if len(oxygen_neighbors) == 2:
                if any(n.GetTotalNumHs() > 0 for n in oxygen_neighbors):
                    hemiacetal_carbons.append(atom)

    if not hemiacetal_carbons:
        return False, "No hemiacetal carbon found"

    # Check if the hemiacetal carbon is part of a monosaccharide or oligosaccharide structure
    for hemiacetal_carbon in hemiacetal_carbons:
        neighbors = hemiacetal_carbon.GetNeighbors()
        if len(neighbors) == 4:  # Check for tetrahedral geometry
            oxygen_neighbors = [n for n in neighbors if n.GetSymbol() == 'O']
            if len(oxygen_neighbors) == 2:
                if (any(n.GetTotalNumHs() == 1 for n in oxygen_neighbors) and
                        any(n.GetTotalNumHs() == 0 for n in oxygen_neighbors)):
                    return True, "Hemiacetal carbon found in a potential glycosyl group"

    return False, "Hemiacetal carbon not part of a glycosyl group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24403',
                          'name': 'glycosyl group',
                          'definition': 'An organic group obtained by removing '
                                        'the hydroxy group from the hemiacetal '
                                        'function of a monosaccharide or '
                                        'monosaccharide derivative and, by '
                                        'extension, of a lower oligosaccharide '
                                        'or oligosaccharide derivative.',
                          'parents': ['CHEBI:33247']},
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 28,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}