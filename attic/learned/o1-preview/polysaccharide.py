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

    # Define a SMARTS pattern for monosaccharide units (pyranose and furanose rings)
    monosaccharide_smarts = "[C;R1]1[C,O][C,O][C,O][C,O][C,O]1"  # 5-membered ring
    monosaccharide_smarts_6 = "[C;R1]1[C,O][C,O][C,O][C,O][C,O][C,O]1"  # 6-membered ring

    monosaccharide_pattern_5 = Chem.MolFromSmarts(monosaccharide_smarts)
    monosaccharide_pattern_6 = Chem.MolFromSmarts(monosaccharide_smarts_6)

    # Find all monosaccharide units
    mono_matches_5 = mol.GetSubstructMatches(monosaccharide_pattern_5)
    mono_matches_6 = mol.GetSubstructMatches(monosaccharide_pattern_6)

    num_monosaccharides = len(mono_matches_5) + len(mono_matches_6)
    if num_monosaccharides == 0:
        return False, "No monosaccharide units found"

    if num_monosaccharides <= 10:
        return False, f"Only {num_monosaccharides} monosaccharide units found, need more than 10"

    # Define a SMARTS pattern for glycosidic linkage (C-O-C between rings)
    glycosidic_bond_smarts = "[C;R][O][C;R]"
    glycosidic_bond_pattern = Chem.MolFromSmarts(glycosidic_bond_smarts)

    # Find glycosidic bonds
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    num_glycosidic_bonds = len(glycosidic_matches)

    if num_glycosidic_bonds < (num_monosaccharides - 1):
        return False, f"Insufficient glycosidic bonds connecting monosaccharide units: found {num_glycosidic_bonds}, expected at least {num_monosaccharides - 1}"

    return True, f"Contains {num_monosaccharides} monosaccharide units connected via glycosidic bonds"