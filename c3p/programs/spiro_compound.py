"""
Classifies: CHEBI:33599 spiro compound
"""
from rdkit import Chem

def is_spiro_compound(smiles: str):
    """
    Determines if a molecule is a spiro compound (a compound having one atom as the only common member of two rings).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a spiro compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Create a dictionary to count the number of rings each atom is part of
    atom_ring_count = {}
    for ring in atom_rings:
        for atom_idx in ring:
            if atom_idx in atom_ring_count:
                atom_ring_count[atom_idx].append(ring)
            else:
                atom_ring_count[atom_idx] = [ring]

    # Check for any atom that is part of exactly two rings
    for atom_idx, rings in atom_ring_count.items():
        if len(rings) == 2:
            # Ensure that the two rings are distinct
            if len(set(rings[0]).intersection(set(rings[1]))) == 1:
                return True, "Atom {} is the only common member of two rings".format(atom_idx)

    return False, "No atom found as the only common member of two rings"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33599',
                          'name': 'spiro compound',
                          'definition': 'A compound having one atom as the '
                                        'only common member of two rings.',
                          'parents': ['CHEBI:33635']},
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
    'num_true_positives': 38,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 2,
    'precision': 1.0,
    'recall': 0.95,
    'f1': 0.9743589743589743,
    'accuracy': None}