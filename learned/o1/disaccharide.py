"""
Classifies: CHEBI:36233 disaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

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

    # Kekulize the molecule to ensure proper valence assignment
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Chem.KekulizeException:
        pass  # Some sugars may contain aromatic rings

    # Define SMARTS patterns for monosaccharide units (pyranose and furanose rings)
    pyranose_smarts = "[C;H1,H2,H3]1[C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][O]1"
    furanose_smarts = "[C;H1,H2,H3]1[C;H1,H2][C;H1,H2][O][C;H1,H2]1"

    pyranose_pattern = Chem.MolFromSmarts(pyranose_smarts)
    furanose_pattern = Chem.MolFromSmarts(furanose_smarts)

    # Identify monosaccharide rings in the molecule
    monosaccharide_rings = []
    for ring in mol.GetRingInfo().AtomRings():
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        submol = Chem.PathToSubmol(mol, ring)
        if submol.HasSubstructMatch(pyranose_pattern) or submol.HasSubstructMatch(furanose_pattern):
            monosaccharide_rings.append(set(ring))

    if len(monosaccharide_rings) != 2:
        return False, f"Found {len(monosaccharide_rings)} monosaccharide units, need exactly 2"

    # Map atom indices to monosaccharide ring indices
    ring_atom_map = {}
    for ring_idx, ring_atoms in enumerate(monosaccharide_rings):
        for atom_idx in ring_atoms:
            ring_atom_map[atom_idx] = ring_idx

    # Find glycosidic bonds connecting the two monosaccharide units
    glycosidic_bonds = []
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        idx1 = atom1.GetIdx()
        idx2 = atom2.GetIdx()
        atomic_nums = (atom1.GetAtomicNum(), atom2.GetAtomicNum())

        # Check for oxygen bridge between rings (C-O-C linkage)
        if 8 in atomic_nums and idx1 in ring_atom_map and idx2 in ring_atom_map:
            ring_idx1 = ring_atom_map[idx1]
            ring_idx2 = ring_atom_map[idx2]
            if ring_idx1 != ring_idx2:
                glycosidic_bonds.append(bond)

    if len(glycosidic_bonds) == 0:
        return False, "No glycosidic bond connecting monosaccharide units found"

    # Ensure that monosaccharide units are connected via a glycosidic bond
    connected_units = set()
    for bond in glycosidic_bonds:
        idx1 = bond.GetBeginAtom().GetIdx()
        idx2 = bond.GetEndAtom().GetIdx()
        ring_idx1 = ring_atom_map[idx1]
        ring_idx2 = ring_atom_map[idx2]
        connected_units.update([ring_idx1, ring_idx2])

    if len(connected_units) != 2:
        return False, "Monosaccharide units are not properly connected via glycosidic bond"

    return True, "Contains two monosaccharide units connected via a glycosidic bond"