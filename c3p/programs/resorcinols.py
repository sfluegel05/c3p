"""
Classifies: CHEBI:33572 resorcinols
"""
from rdkit import Chem

def is_resorcinols(smiles: str):
    """
    Determines if a molecule is a resorcinol (benzenediol with two hydroxy groups meta to one another).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a resorcinol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least one 6-membered ring
    if not any(len(ring) == 6 for ring in rings.AtomRings()):
        return False, "No 6-membered rings found"

    # Find all 6-membered rings with exactly two hydroxy groups
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            hydroxy_positions = [i for i, atom in enumerate(atoms) if atom.GetSymbol() == 'C' and any(neigh.GetSymbol() == 'O' and neigh.GetTotalNumHs() == 1 for neigh in atom.GetNeighbors())]
            if len(hydroxy_positions) == 2:
                # Check if the hydroxy groups are meta to each other
                pos1, pos2 = hydroxy_positions
                if abs(pos1 - pos2) == 2 or abs(pos1 - pos2) == 4:
                    return True, "Molecule is a resorcinol"
    
    return False, "Molecule is not a resorcinol"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33572',
                          'name': 'resorcinols',
                          'definition': 'Any benzenediol in which the two '
                                        'hydroxy groups are meta to one '
                                        'another.',
                          'parents': ['CHEBI:33570']},
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
    'num_true_positives': 25,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 4,
    'precision': 0.9259259259259259,
    'recall': 0.8620689655172413,
    'f1': 0.8928571428571429,
    'accuracy': None}