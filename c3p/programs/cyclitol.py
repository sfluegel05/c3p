"""
Classifies: CHEBI:23451 cyclitol
"""
from rdkit import Chem

def is_cyclitol(smiles: str):
    """
    Determines if a molecule is a cyclitol (a polyol consisting of a cycloalkane containing at least three hydroxy groups, each attached to a different ring carbon atom).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclitol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for cycloalkane rings
    cycloalkane_rings = []
    for ring in rings.AtomRings():
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if all(atom.GetSymbol() == 'C' for atom in atoms):
            cycloalkane_rings.append(ring)

    if not cycloalkane_rings:
        return False, "No cycloalkane rings found"

    # Check for at least three hydroxy groups attached to different ring carbon atoms
    for ring in cycloalkane_rings:
        hydroxy_count = 0
        hydroxy_positions = set()
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O' and neighbor.GetTotalDegree() == 1:
                    hydroxy_count += 1
                    hydroxy_positions.add(atom_idx)
        if hydroxy_count >= 3 and len(hydroxy_positions) >= 3:
            return True, "Molecule is a cyclitol"

    return False, "Molecule does not meet the cyclitol criteria"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23451',
                          'name': 'cyclitol',
                          'definition': 'A polyol consisting of a cycloalkane '
                                        'containing at least three hydroxy '
                                        'groups, each attached to a different '
                                        'ring carbon atom.',
                          'parents': ['CHEBI:26191']},
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
    'num_false_negatives': 31,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}