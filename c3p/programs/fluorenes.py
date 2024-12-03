"""
Classifies: CHEBI:24059 fluorenes
"""
from rdkit import Chem

def is_fluorenes(smiles: str):
    """
    Determines if a molecule is a fluorene.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fluorene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()
    atom_rings = rings.AtomRings()

    # Check for at least three rings
    if len(atom_rings) < 3:
        return False, "Less than three rings found"

    # Find all 6-membered rings and 5-membered rings
    six_membered_rings = [ring for ring in atom_rings if len(ring) == 6]
    five_membered_rings = [ring for ring in atom_rings if len(ring) == 5]

    if len(six_membered_rings) < 2:
        return False, "Less than two 6-membered rings found"

    if len(five_membered_rings) < 1:
        return False, "No 5-membered rings found"

    # Check for ortho-fusion between the benzene rings and the cyclopentane ring
    for five_ring in five_membered_rings:
        five_ring_set = set(five_ring)
        ortho_fused = 0
        for six_ring in six_membered_rings:
            six_ring_set = set(six_ring)
            if len(five_ring_set.intersection(six_ring_set)) >= 2:
                ortho_fused += 1
        if ortho_fused >= 2:
            return True, "Fluorene structure found"

    return False, "No ortho-fused structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24059',
                          'name': 'fluorenes',
                          'definition': 'An ortho-fused polycyclic arene in '
                                        'which the skeleton is composed of two '
                                        'benzene rings ortho-fused to '
                                        'cyclopentane.',
                          'parents': ['CHEBI:38032']},
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
    'num_true_positives': 10,
    'num_false_positives': 1,
    'num_true_negatives': 9,
    'num_false_negatives': 0,
    'precision': 0.9090909090909091,
    'recall': 1.0,
    'f1': 0.9523809523809523,
    'accuracy': None}