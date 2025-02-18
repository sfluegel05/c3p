"""
Classifies: CHEBI:50699 oligosaccharide
"""
"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles

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

    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    
    # Find all 5 or 6-membered oxygen-containing rings (potential monosaccharide units)
    sugar_rings = []
    for ring in ring_info.AtomRings():
        ring_size = len(ring)
        if ring_size not in (5, 6):
            continue
        # Check if ring contains an oxygen atom
        if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring):
            sugar_rings.append(set(ring))

    if len(sugar_rings) < 2:
        return False, f"Found {len(sugar_rings)} sugar rings, need at least 2"

    # Count glycosidic bonds (bridging oxygens between different rings)
    glycosidic_bonds = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:
            neighbors = atom.GetNeighbors()
            if len(neighbors) != 2:
                continue  # Not part of a linkage

            # Check if neighbors belong to different sugar rings
            ring_membership = []
            for neighbor in neighbors:
                for i, s_ring in enumerate(sugar_rings):
                    if neighbor.GetIdx() in s_ring:
                        ring_membership.append(i)
                        break  # Each atom can only belong to one sugar ring

            if len(ring_membership) == 2 and ring_membership[0] != ring_membership[1]:
                glycosidic_bonds += 1

    required_bonds = len(sugar_rings) - 1
    if glycosidic_bonds < required_bonds:
        return False, f"Found {glycosidic_bonds} glycosidic bonds, need at least {required_bonds}"

    return True, f"Contains {len(sugar_rings)} sugar units with {glycosidic_bonds} glycosidic bonds"