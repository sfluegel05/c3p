"""
Classifies: CHEBI:36233 disaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is a compound in which two monosaccharides are joined by a glycosidic bond.

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

    # Define SMARTS patterns for pyranose and furanose rings
    # These patterns match common sugar rings, accounting for possible substitutions
    pyranose_smarts = "[#6]-1-[#6]-[#6]-[#6]-[#6]-[#8]-1"
    furanose_smarts = "[#6]-1-[#6]-[#6]-[#6]-[#8]-1"

    pyranose = Chem.MolFromSmarts(pyranose_smarts)
    furanose = Chem.MolFromSmarts(furanose_smarts)

    # Identify sugar rings using the SMARTS patterns
    sugar_rings = []
    for ring in mol.GetRingInfo().AtomRings():
        submol = Chem.PathToSubmol(mol, ring)
        if submol.HasSubstructMatch(pyranose) or submol.HasSubstructMatch(furanose):
            sugar_rings.append(set(ring))

    if len(sugar_rings) < 2:
        return False, f"Found {len(sugar_rings)} sugar rings, need at least 2"

    # Map atom indices to sugar ring indices
    ring_atom_map = {}
    for ring_idx, ring_atoms in enumerate(sugar_rings):
        for atom_idx in ring_atoms:
            ring_atom_map[atom_idx] = ring_idx

    # Find glycosidic bonds connecting two sugar rings
    glycosidic_bonds = []
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        idx1 = atom1.GetIdx()
        idx2 = atom2.GetIdx()

        # Check if bond connects two sugar rings
        in_ring1 = idx1 in ring_atom_map
        in_ring2 = idx2 in ring_atom_map
        if in_ring1 and in_ring2:
            ring_idx1 = ring_atom_map[idx1]
            ring_idx2 = ring_atom_map[idx2]
            if ring_idx1 != ring_idx2:
                # Bond connects two different sugar rings
                glycosidic_bonds.append((idx1, idx2))
        else:
            # Check for oxygen atoms connecting two sugar rings
            if atom1.GetAtomicNum() == 8 and idx2 in ring_atom_map:
                for neighbor in atom1.GetNeighbors():
                    n_idx = neighbor.GetIdx()
                    if n_idx != idx2 and n_idx in ring_atom_map and ring_atom_map[n_idx] != ring_atom_map[idx2]:
                        glycosidic_bonds.append((n_idx, idx1, idx2))
            elif atom2.GetAtomicNum() == 8 and idx1 in ring_atom_map:
                for neighbor in atom2.GetNeighbors():
                    n_idx = neighbor.GetIdx()
                    if n_idx != idx1 and n_idx in ring_atom_map and ring_atom_map[n_idx] != ring_atom_map[idx1]:
                        glycosidic_bonds.append((n_idx, idx2, idx1))

    if len(glycosidic_bonds) == 0:
        return False, "No glycosidic bond connecting sugar rings found"

    # Count the number of monosaccharide units
    monosaccharide_units = set(ring_atom_map.values())
    if len(monosaccharide_units) != 2:
        return False, f"Found {len(monosaccharide_units)} monosaccharide units, need exactly 2"

    # Ensure that the monosaccharide units are connected via glycosidic bonds
    connected_units = set()
    for bond in glycosidic_bonds:
        # bond can be (idx1, idx2) or (n_idx, o_idx, idx)
        if len(bond) == 2:
            idx1, idx2 = bond
            ring_idx1 = ring_atom_map[idx1]
            ring_idx2 = ring_atom_map[idx2]
        else:
            idx1, _, idx2 = bond
            ring_idx1 = ring_atom_map[idx1]
            ring_idx2 = ring_atom_map[idx2]
        connected_units.update([ring_idx1, ring_idx2])

    if len(connected_units) != 2:
        return False, "Monosaccharide units are not properly connected via glycosidic bonds"

    return True, "Contains two monosaccharide units connected via a glycosidic bond"