"""
Classifies: CHEBI:51681 dimethoxybenzene
"""
from rdkit import Chem

def is_dimethoxybenzene(smiles: str):
    """
    Determines if a molecule is a dimethoxybenzene (benzene substituted with two methoxy groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dimethoxybenzene, False otherwise
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

    # Find all aromatic 6-membered rings
    aromatic_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() for atom in atoms):
                aromatic_rings.append(ring)

    if not aromatic_rings:
        return False, "No aromatic 6-membered rings found"

    # Check for methoxy groups
    methoxy_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and len(atom.GetNeighbors()) == 1:
            neighbor = atom.GetNeighbors()[0]
            if neighbor.GetSymbol() == 'C' and len(neighbor.GetNeighbors()) == 3:
                if any(neigh.GetSymbol() == 'C' and neigh.GetIsAromatic() for neigh in neighbor.GetNeighbors()):
                    methoxy_count += 1

    if methoxy_count < 2:
        return False, "Less than two methoxy groups found"

    return True, "Dimethoxybenzene identified"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51681',
                          'name': 'dimethoxybenzene',
                          'definition': 'Any methoxybenzene that consists of a '
                                        'benzene skeleton substituted with two '
                                        'methoxy groups and its derivatives.',
                          'parents': ['CHEBI:51683']},
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
    'num_true_positives': 9,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 24,
    'precision': 0.8181818181818182,
    'recall': 0.2727272727272727,
    'f1': 0.4090909090909091,
    'accuracy': None}