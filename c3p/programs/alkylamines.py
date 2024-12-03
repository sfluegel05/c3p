"""
Classifies: CHEBI:22331 alkylamines
"""
from rdkit import Chem

def is_alkylamines(smiles: str):
    """
    Determines if a molecule is an alkylamine (amine derived from ammonia by replacing one, two, or three hydrogen atoms by alkyl groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkylamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of nitrogen atoms
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'N']
    if not nitrogen_atoms:
        return False, "No nitrogen atoms found"

    # Check if nitrogen atoms are connected to alkyl groups
    for nitrogen in nitrogen_atoms:
        connected_atoms = nitrogen.GetNeighbors()
        non_hydrogen_neighbors = [atom for atom in connected_atoms if atom.GetSymbol() != 'H']
        if len(non_hydrogen_neighbors) > 3:
            continue  # Not an amine if nitrogen has more than 3 non-hydrogen neighbors

        alkyl_groups = 0
        for neighbor in non_hydrogen_neighbors:
            if neighbor.GetSymbol() == 'C' and all(neigh.GetSymbol() in ['C', 'H'] for neigh in neighbor.GetNeighbors() if neigh.GetIdx() != nitrogen.GetIdx()):
                alkyl_groups += 1

        if alkyl_groups >= 1:
            return True, "Molecule is an alkylamine"

    return False, "No alkylamine structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22331',
                          'name': 'alkylamines',
                          'definition': 'Any amine formally derived from '
                                        'ammonia by replacing one, two or '
                                        'three hydrogen atoms by alkyl groups.',
                          'parents': ['CHEBI:32952']},
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
    'num_true_positives': 38,
    'num_false_positives': 18,
    'num_true_negatives': 2,
    'num_false_negatives': 0,
    'precision': 0.6785714285714286,
    'recall': 1.0,
    'f1': 0.8085106382978724,
    'accuracy': None}