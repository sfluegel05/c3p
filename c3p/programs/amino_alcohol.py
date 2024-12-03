"""
Classifies: CHEBI:22478 amino alcohol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def is_amino_alcohol(smiles: str):
    """
    Determines if a molecule is an amino alcohol (an alcohol containing an amino functional group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of at least one hydroxyl group (OH)
    hydroxyl_groups = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and any(neigh.GetSymbol() == 'C' for neigh in atom.GetNeighbors())]
    if not hydroxyl_groups:
        return False, "No hydroxyl groups found"

    # Check for the presence of at least one amino group (NH2, NH, or N)
    amino_groups = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'N']
    if not amino_groups:
        return False, "No amino groups found"

    # Ensure there is a carbon atom bonded to both a hydroxyl and an amino group
    for hydroxyl in hydroxyl_groups:
        for neighbor in hydroxyl.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                for neigh in neighbor.GetNeighbors():
                    if neigh.GetSymbol() == 'N':
                        return True, "Molecule is an amino alcohol"
    
    return False, "No carbon found bonded to both hydroxyl and amino groups"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22478',
                          'name': 'amino alcohol',
                          'definition': 'An alcohol containing an amino '
                                        'functional group in addition to the '
                                        'alcohol-defining hydroxy group.',
                          'parents': ['CHEBI:30879', 'CHEBI:50047']},
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
    'num_true_positives': 12,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 21,
    'precision': 0.8571428571428571,
    'recall': 0.36363636363636365,
    'f1': 0.5106382978723404,
    'accuracy': None}