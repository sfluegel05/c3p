"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: disaccharide
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is composed of two monosaccharide units joined by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    ring_atoms = ring_info.AtomRings()
    sugar_rings = []
    sugar_ring_indices = []  # Keep track of ring indices

    # Identify sugar rings
    for idx, ring in enumerate(ring_atoms):
        # Consider rings of size 5 or 6
        if len(ring) == 5 or len(ring) == 6:
            is_sugar_ring = True
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                # Check if atom is carbon or oxygen
                if atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 8:
                    is_sugar_ring = False
                    break
                # Check that bonds are single
                for bond in atom.GetBonds():
                    if bond.GetBondTypeAsDouble() != 1.0:
                        is_sugar_ring = False
                        break
                if not is_sugar_ring:
                    break
            if is_sugar_ring:
                sugar_rings.append(set(ring))
                sugar_ring_indices.append(idx)

    # Ensure there are exactly two sugar rings
    if len(sugar_rings) != 2:
        return False, f"Found {len(sugar_rings)} sugar rings, need exactly 2"

    # Create atom index to sugar ring index mapping
    atom_to_ring = {}
    for ring_idx, ring in enumerate(sugar_rings):
        for atom_idx in ring:
            atom_to_ring[atom_idx] = ring_idx

    # Identify glycosidic bonds between sugar rings
    glycosidic_bonds = []
    for atom in mol.GetAtoms():
        # Non-ring oxygen atom of degree 2 (connected to two atoms)
        if atom.GetAtomicNum() == 8 and not atom.IsInRing() and atom.GetDegree() == 2:
            neighbors = atom.GetNeighbors()
            if len(neighbors) != 2:
                continue
            # Get ring indices of neighboring atoms
            neighbor_rings = []
            for neighbor in neighbors:
                neighbor_idx = neighbor.GetIdx()
                ring_idx = atom_to_ring.get(neighbor_idx, None)
                neighbor_rings.append(ring_idx)
            # Check if connected to atoms from two different sugar rings
            if neighbor_rings[0] is not None and neighbor_rings[1] is not None:
                if neighbor_rings[0] != neighbor_rings[1]:
                    glycosidic_bonds.append((atom.GetIdx(), neighbors[0].GetIdx(), neighbors[1].GetIdx()))

    # Ensure there is exactly one glycosidic bond connecting the two sugar rings
    if len(glycosidic_bonds) != 1:
        return False, f"Found {len(glycosidic_bonds)} glycosidic bonds connecting sugar rings, need exactly 1"

    # Check for additional sugar rings (to exclude oligosaccharides)
    total_rings = len(sugar_rings)
    if total_rings != 2:
        return False, f"Found {total_rings} sugar rings, indicates more than two sugar units"

    return True, "Contains two monosaccharide units connected by a glycosidic bond"

__metadata__ = {   'chemical_class': {   'id': None,
                          'name': 'disaccharide',
                          'definition': 'A compound in which two monosaccharides are joined by a glycosidic bond.',
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
    'attempt': 3,
    'success': None,
    'best': None,
    'error': '',
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None}