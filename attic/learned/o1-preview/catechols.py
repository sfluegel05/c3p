"""
Classifies: CHEBI:33566 catechols
"""
"""
Classifies: catechols
"""

from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule contains a catechol moiety based on its SMILES string.
    A catechol moiety is an o-diphenol component: an aromatic ring with two adjacent hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a catechol moiety, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    # Iterate over each ring
    for ring in rings:
        # Check if ring is aromatic
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            ring_size = len(ring)
            # Iterate over atoms in the ring
            for i in range(ring_size):
                idx1 = ring[i]
                idx2 = ring[(i+1)%ring_size]  # Next atom in the ring (with wrap-around)
                atom1 = mol.GetAtomWithIdx(idx1)
                atom2 = mol.GetAtomWithIdx(idx2)
                # Check if both atoms are aromatic carbons
                if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                    # Check if atom1 has -OH group
                    has_OH1 = False
                    for neighbor in atom1.GetNeighbors():
                        if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                            has_OH1 = True
                            break
                    # Check if atom2 has -OH group
                    has_OH2 = False
                    for neighbor in atom2.GetNeighbors():
                        if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                            has_OH2 = True
                            break
                    # If both have -OH groups, it's a catechol moiety
                    if has_OH1 and has_OH2:
                        return True, "Contains o-diphenol component (catechol moiety)"
    # If no catechol moiety found
    return False, "No o-diphenol component found"

__metadata__ = {   'chemical_class': {   'name': 'catechols',
                              'definition': 'Any compound containing an o-diphenol component.'},
        'config': {   'llm_model_name': 'lbl/claude-sonnet',
                      'f1_threshold': 0.8,
                      'max_attempts': 5,
                      'max_positive_instances': None,
                      'max_positive_to_test': None,
                      'max_negative_to_test': None,
                      'max_positive_in_prompt': 50,
                      'max_negative_in_prompt': 20,
                      'max_instances_in_prompt': 100,
                      'test_proportion': 0.1},
        'message': None,
        'attempt': 1,
        'success': True,
        'best': True,
        'error': '',
        'stdout': None}