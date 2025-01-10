"""
Classifies: CHEBI:72600 spiroketal
"""
"""
Classifies: CHEBI:35567 spiroketal
"""
from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal is a cyclic ketal in which the ketal carbon is the only common atom of two rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a spiroketal, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Find spiro atoms in the molecule
    spiro_atoms = Chem.FindSpiroAtoms(mol)

    # Iterate over spiro atoms to check for spiroketal centers
    for spiro_idx in spiro_atoms:
        atom = mol.GetAtomWithIdx(spiro_idx)
        if atom.GetAtomicNum() != 6:
            continue  # Skip if not a carbon atom

        if atom.GetDegree() != 4:
            continue  # Spiroketal carbon should have degree 4

        neighbors = atom.GetNeighbors()
        oxy_count = 0
        oxy_ring_count = 0
        for nbr in neighbors:
            if nbr.GetAtomicNum() == 8:
                oxy_count += 1
                if nbr.IsInRing():
                    oxy_ring_count += 1
        if oxy_count != 2:
            continue  # Spiroketal carbon should be connected to exactly two oxygens
        if oxy_ring_count != 2:
            continue  # Both oxygens should be part of rings

        # Check that the spiro carbon is common to exactly two rings that share only this atom
        rings_with_atom = [ring for ring in atom_rings if spiro_idx in ring]
        if len(rings_with_atom) != 2:
            continue  # Spiroketal carbon should be in exactly two rings

        # Check that the two rings share only the spiro carbon atom
        shared_atoms = set(rings_with_atom[0]) & set(rings_with_atom[1])
        if len(shared_atoms) != 1:
            continue  # Rings share more than the spiro carbon

        return True, f"Spiroketal detected at atom index {spiro_idx}"

    return False, "No spiroketal center found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35567',
                                  'name': 'spiroketal',
                                  'definition': 'A cyclic ketal in which the ketal '
                                                'carbon is the only common atom of '
                                                'two rings.',
                                  'parents': []},
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
        'attempt': 2,
        'success': True,
        'best': True,
        'error': '',
        'stdout': None
    }