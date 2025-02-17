"""
Classifies: CHEBI:50126 tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide is defined as an oligosaccharide comprising four monosaccharide units.
    In many cases, each monosaccharide unit appears as a cyclic structure (typically a furanose or pyranose)
    that contains one ring oxygen with the remaining atoms being carbons.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a tetrasaccharide, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES input
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Define a helper function to test if a ring is a typical monosaccharide ring.
    def is_monosaccharide_ring(atom_indices):
        """
        Returns True if the ring defined by atom_indices is a candidate monosaccharide ring.
        We assume that a furanose ring (5-membered) should contain 1 oxygen and 4 carbons,
        and a pyranose ring (6-membered) should contain 1 oxygen and 5 carbons.
        """
        ring_size = len(atom_indices)
        if ring_size not in [5, 6]:
            return False
        atoms = [mol.GetAtomWithIdx(idx) for idx in atom_indices]
        o_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
        c_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
        if ring_size == 5:
            return (o_count == 1 and c_count == 4)
        if ring_size == 6:
            return (o_count == 1 and c_count == 5)
        return False

    # Count valid monosaccharide rings
    valid_ring_count = 0
    for ring in atom_rings:
        if is_monosaccharide_ring(ring):
            valid_ring_count += 1

    # Classify as a tetrasaccharide if exactly 4 monosaccharide rings are found.
    if valid_ring_count == 4:
        return True, "Contains exactly four monosaccharide rings typical for tetrasaccharides."
    else:
        return False, f"Found {valid_ring_count} monosaccharide-like rings; expected 4 for a tetrasaccharide."