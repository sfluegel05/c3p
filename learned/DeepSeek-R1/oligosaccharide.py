"""
Classifies: CHEBI:50699 oligosaccharide
"""
"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, GetRingInfo

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    An oligosaccharide consists of multiple monosaccharide units joined by glycosidic linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligosaccharide, False otherwise
        str: Reason for classification
    """
    mol = MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all 5 or 6-membered rings containing oxygen (potential monosaccharide units)
    ring_info = GetRingInfo(mol)
    sugar_rings = []
    for ring in ring_info.AtomRings():
        ring_size = len(ring)
        if ring_size not in (5, 6):
            continue
        has_oxygen = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring)
        if has_oxygen:
            sugar_rings.append(set(ring))

    if len(sugar_rings) < 2:
        return False, f"Found {len(sugar_rings)} sugar rings, need at least 2"

    # Count bridging oxygen atoms between different sugar rings
    bridging_os = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:
            neighbors = atom.GetNeighbors()
            if len(neighbors) != 2:
                continue  # Not part of a glycosidic bond

            # Check if neighbors are in different sugar rings
            neighbor_rings = []
            for neighbor in neighbors:
                for i, s_ring in enumerate(sugar_rings):
                    if neighbor.GetIdx() in s_ring:
                        neighbor_rings.append(i)
                        break  # Assume each neighbor is in only one sugar ring

            if len(neighbor_rings) == 2 and neighbor_rings[0] != neighbor_rings[1]:
                bridging_os += 1

    required_bonds = len(sugar_rings) - 1
    if bridging_os < required_bonds:
        return False, f"Found {bridging_os} glycosidic bonds, need at least {required_bonds}"

    return True, f"Contains {len(sugar_rings)} monosaccharide units with {bridging_os} glycosidic bonds"