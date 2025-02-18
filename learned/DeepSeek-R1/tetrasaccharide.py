"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: CHEBI:166592 tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide is an oligosaccharide comprising four monomeric monosaccharide units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"

    # Find all rings that are potential sugar rings (5 or 6-membered with oxygen)
    ri = mol.GetRingInfo()
    sugar_rings = []
    for ring in ri.AtomRings():
        ring_size = len(ring)
        if ring_size not in [5, 6]:
            continue
        # Check if the ring contains an oxygen atom
        has_oxygen = any(mol.GetAtomWithIdx(a).GetAtomicNum() == 8 for a in ring)
        if has_oxygen:
            sugar_rings.append(ring)

    if len(sugar_rings) != 4:
        return False, f"Found {len(sugar_rings)} sugar rings, expected 4"

    # Find glycosidic bonds (ethers connecting different sugar rings)
    glycosidic_bonds = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetDegree() == 2:
            neighbors = atom.GetNeighbors()
            if len(neighbors) != 2:
                continue
            a1, a2 = neighbors
            # Determine which sugar rings the neighbors belong to
            ring_indices = []
            for a in [a1, a2]:
                for i, r in enumerate(sugar_rings):
                    if a.GetIdx() in r:
                        ring_indices.append(i)
                        break
            if len(ring_indices) == 2 and ring_indices[0] != ring_indices[1]:
                glycosidic_bonds += 1

    if glycosidic_bonds < 3:
        return False, f"Found {glycosidic_bonds} glycosidic bonds, need at least 3"

    # Additional check for molecular formula (optional)
    # Calculate the number of oxygen atoms (each sugar has multiple O's, but this varies)
    # Not reliable, but could add if needed.

    return True, "Four sugar rings connected by three glycosidic bonds"