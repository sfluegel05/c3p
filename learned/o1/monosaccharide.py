"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: monosaccharide
Determines if a molecule is a monosaccharide based on its SMILES string.
"""

from rdkit import Chem

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    A monosaccharide is a polyhydroxy aldehyde or ketone with three or more carbon atoms,
    without glycosidic connections to other saccharide units. It includes cyclic forms
    (furanose and pyranose), deoxy sugars, uronic acids, and their derivatives, provided
    that the parent compound has a (potential) carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check number of carbon atoms (must be 3 or more)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return False, f"Contains {c_count} carbon atoms, which is less than 3"

    # Check allowed elements (C, H, O, N, P)
    allowed_atomic_nums = {1, 6, 7, 8, 15}  # H, C, N, O, P
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains element with atomic number {atom.GetAtomicNum()}, which is not typical in monosaccharides"

    # Check for glycosidic bonds (acetal linkages at anomeric carbon)
    # Anomeric carbon in glycosidic bond: a carbon connected to two oxygens via single bonds (C-O-C-O)
    glycosidic_pattern = Chem.MolFromSmarts("[C;R0;X4](O)(O)[#6]")
    if mol.HasSubstructMatch(glycosidic_pattern):
        return False, "Contains glycosidic bonds indicating connection to other saccharide units"

    # Detect furanose ring (5-membered ring with 4 carbons and 1 oxygen)
    furanose_pattern = Chem.MolFromSmarts("C1OC(CO)C1O")
    has_furanose = mol.HasSubstructMatch(furanose_pattern)

    # Detect pyranose ring (6-membered ring with 5 carbons and 1 oxygen)
    pyranose_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1O")
    has_pyranose = mol.HasSubstructMatch(pyranose_pattern)

    # Detect open-chain aldose (polyhydroxy aldehyde)
    aldose_pattern = Chem.MolFromSmarts("[CH](=O)[CH2][CH](O)[CH](O)[CH2]O")
    has_aldose = mol.HasSubstructMatch(aldose_pattern)

    # Detect open-chain ketose (polyhydroxy ketone)
    ketose_pattern = Chem.MolFromSmarts("[CH2]O[CH](O)[CH](O)[C](=O)[CH2]O")
    has_ketose = mol.HasSubstructMatch(ketose_pattern)

    # Check for at least two hydroxyl groups (-OH)
    hydroxyl_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() > 0]
    if len(hydroxyl_atoms) < 2:
        return False, f"Contains {len(hydroxyl_atoms)} hydroxyl groups, which is less than 2"

    # Check for carbonyl group (aldehyde or ketone)
    carbonyl_pattern = Chem.MolFromSmarts("C=O")
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)

    # Check for hemiacetal or hemiketal (cyclic forms)
    hemiacetal_pattern = Chem.MolFromSmarts("[C;R][O][C;R]")
    has_hemiacetal = mol.HasSubstructMatch(hemiacetal_pattern)

    # Monosaccharide should have either cyclic form or open-chain form with carbonyl
    if not (has_furanose or has_pyranose or has_aldose or has_ketose or has_hemiacetal or has_carbonyl):
        return False, "Does not contain structural features typical of monosaccharides"

    # Exclude molecules with multiple sugar rings (polysaccharides)
    ri = mol.GetRingInfo()
    ring_atoms = ri.AtomRings()
    sugar_rings = 0
    for ring in ring_atoms:
        ring_size = len(ring)
        o_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if (ring_size == 5 or ring_size == 6) and o_in_ring == 1:
            sugar_rings += 1
    if sugar_rings > 1:
        return False, f"Contains {sugar_rings} sugar rings, which is more than 1"

    # Check molecule size (exclude large molecules)
    heavy_atom_count = mol.GetNumHeavyAtoms()
    if heavy_atom_count > 30:
        return False, f"Contains {heavy_atom_count} heavy atoms, which is too large for a monosaccharide"

    return True, "Molecule satisfies criteria for a monosaccharide"