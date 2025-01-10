"""
Classifies: CHEBI:18154 polysaccharide
"""
"""
Classifies: polysaccharide
"""
from rdkit import Chem

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide is defined as a biomacromolecule consisting of large numbers
    of monosaccharide residues linked glycosidically, typically containing more
    than ten monosaccharide units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Identify monosaccharide rings (5- or 6-membered rings with 1 oxygen and rest carbons)
    num_monosaccharides = 0
    monosaccharide_rings = []

    for ring in atom_rings:
        ring_size = len(ring)
        if ring_size == 5 or ring_size == 6:
            num_oxygen = 0
            num_non_carbon = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                atomic_num = atom.GetAtomicNum()
                if atomic_num == 8:
                    num_oxygen += 1
                if atomic_num != 6:
                    num_non_carbon += 1
            if num_oxygen == 1 and num_non_carbon == 1:
                # Found a monosaccharide ring
                num_monosaccharides += 1
                monosaccharide_rings.append(set(ring))

    if num_monosaccharides == 0:
        return False, "No monosaccharide units found"

    if num_monosaccharides <= 10:
        return False, f"Only {num_monosaccharides} monosaccharide units found, need more than 10"

    # Check for glycosidic bonds between monosaccharide units
    # Glycosidic bonds are ether linkages (C-O-C) connecting different rings
    glycosidic_bonds = 0

    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        atomic_num1 = atom1.GetAtomicNum()
        atomic_num2 = atom2.GetAtomicNum()
        idx1 = atom1.GetIdx()
        idx2 = atom2.GetIdx()

        # Look for C-O single bonds
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            if (atomic_num1 == 6 and atomic_num2 == 8) or (atomic_num1 == 8 and atomic_num2 == 6):
                # Check if this bond connects two different monosaccharide rings
                rings1 = [ring for ring in monosaccharide_rings if idx1 in ring]
                rings2 = [ring for ring in monosaccharide_rings if idx2 in ring]
                if rings1 and rings2 and rings1 != rings2:
                    # Ether linkage between two monosaccharide units
                    glycosidic_bonds += 1

    expected_bonds = num_monosaccharides - 1
    if glycosidic_bonds < expected_bonds:
        return False, f"Insufficient glycosidic bonds connecting monosaccharide units: found {glycosidic_bonds}, expected at least {expected_bonds}"

    return True, f"Contains {num_monosaccharides} monosaccharide units connected via glycosidic bonds"