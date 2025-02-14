"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: monosaccharide
Determines if a molecule is a monosaccharide based on its SMILES string.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    A monosaccharide is a polyhydroxy aldehyde or ketone with three or more carbon atoms,
    without glycosidic connections to other saccharide units. It includes cyclic hemiacetals,
    hemiketals (furanose and pyranose forms), and their derivatives like deoxy sugars.

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

    # Check for aldehyde group (terminal carbonyl group)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)

    # Check for ketone group (internal carbonyl group)
    ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    has_ketone = mol.HasSubstructMatch(ketone_pattern)

    # Check for hemiacetal or hemiketal functionality
    # Look for carbon bonded to -OH and -O- (hemiacetal center)
    hemiacetal_pattern = Chem.MolFromSmarts("[CX4H1]([OX2H])[OX2]")
    has_hemiacetal = mol.HasSubstructMatch(hemiacetal_pattern)

    # Monosaccharides should have at least one of the above functionalities
    if not (has_aldehyde or has_ketone or has_hemiacetal):
        return False, "Does not contain an aldehyde, ketone, or hemiacetal/hemiketal group"

    # Check for multiple hydroxyl groups (-OH) attached to sp3 carbons
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4;H1,H2][OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, f"Contains {len(hydroxyl_matches)} hydroxyl groups attached to carbons, which is less than 1"

    # Ensure molecule is a single saccharide unit (no glycosidic bonds)
    # Glycosidic bonds are acetal linkages at anomeric carbon connecting two saccharides
    # Anomeric carbon is acetal (connected to two oxygens) in glycosidic bonds
    glycosidic_pattern = Chem.MolFromSmarts("[CX4H0](O)(O)[#6]")
    has_glycosidic_bond = mol.HasSubstructMatch(glycosidic_pattern)
    if has_glycosidic_bond:
        return False, "Contains glycosidic bonds indicating connection to other saccharide units"

    # Check for carboxylic acids or esters, which are common in uronic acids
    # Include these derivatives in monosaccharides
    # Uronic acids have a carboxylic acid group at the terminal carbon
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    has_carboxylic_acid = mol.HasSubstructMatch(carboxylic_acid_pattern)

    # Check that the molecule does not have peptide bonds, long alkyl chains, or other non-saccharide components
    # Count total number of non-hydrogen atoms (Exclude very large molecules)
    heavy_atom_count = mol.GetNumHeavyAtoms()
    if heavy_atom_count > 50:
        return False, f"Contains {heavy_atom_count} heavy atoms, which is too large for a monosaccharide"

    # Check that the molecule is primarily composed of C, H, O, and N (allowing for deoxy sugars and amino sugars)
    allowed_atomic_nums = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53}  # H, C, N, O, F, P, S, Cl, Br, I
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains element with atomic number {atom.GetAtomicNum()}, which is not typical in monosaccharides"

    return True, "Molecule satisfies criteria for a monosaccharide"