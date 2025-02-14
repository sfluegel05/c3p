"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: CHEBI:28053 tetrasaccharide
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
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    ring_atoms = ring_info.AtomRings()

    monosaccharide_rings = []
    for ring in ring_atoms:
        if len(ring) in [5, 6]:  # furanose or pyranose
            atom_types = [mol.GetAtomWithIdx(idx).GetAtomicNum() for idx in ring]
            o_in_ring = sum(1 for num in atom_types if num == 8)
            if o_in_ring == 1:  # Ring contains exactly one oxygen atom
                monosaccharide_rings.append(ring)

    # Remove duplicate rings (possible in fused ring systems)
    monosaccharide_rings = [list(ring) for ring in set(tuple(r) for r in monosaccharide_rings)]

    # Check if there are exactly four monosaccharide units
    if len(monosaccharide_rings) != 4:
        return False, f"Found {len(monosaccharide_rings)} monosaccharide units, need exactly 4"

    # Check for glycosidic linkages between the monosaccharide units
    # Glycosidic bonds are typically C-O-C connections between rings
    glycosidic_bonds = 0
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        # Check if bond connects two rings via an oxygen atom
        if (begin_atom.GetAtomicNum() == 8 or end_atom.GetAtomicNum() == 8):
            begin_ring = any(begin_atom.GetIdx() in ring for ring in monosaccharide_rings)
            end_ring = any(end_atom.GetIdx() in ring for ring in monosaccharide_rings)
            if begin_ring and end_ring:
                glycosidic_bonds += 1

    # For a tetrasaccharide, there should be at least 3 glycosidic bonds connecting the 4 units
    if glycosidic_bonds < 3:
        return False, f"Found {glycosidic_bonds} glycosidic bonds, need at least 3"

    return True, "Contains four monosaccharide units linked via glycosidic bonds"