"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: CHEBI:28053 tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Define SMARTS patterns for monosaccharide rings (both furanose and pyranose forms)
    monosaccharide_patterns = [
        # Pyranose form (6-membered ring with one oxygen)
        Chem.MolFromSmarts('C1[C,H][C,H][C,H][C,H]O1'),
        # Furanose form (5-membered ring with one oxygen)
        Chem.MolFromSmarts('C1[C,H][C,H][C,H]O1'),
    ]

    # Find monosaccharide units
    monosaccharide_matches = []
    for pattern in monosaccharide_patterns:
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            # Check if the ring atoms form a ring
            if mol.GetRingInfo().IsAtomInRingOfSize(match[0], len(match)):
                # To avoid duplicates
                if set(match) not in [set(m) for m in monosaccharide_matches]:
                    monosaccharide_matches.append(match)

    # Count unique monosaccharide units
    num_monosaccharides = len(monosaccharide_matches)

    if num_monosaccharides != 4:
        return False, f"Found {num_monosaccharides} monosaccharide units, need exactly 4"

    # Check for glycosidic linkages
    # Glycosidic bond: an ether bond between two monosaccharides
    glycosidic_bonds = 0
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        # Check if the bond is between two oxygen-linked carbons (possible glycosidic bond)
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                begin_oxygens = [nbr for nbr in begin_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
                end_oxygens = [nbr for nbr in end_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
                if begin_oxygens and end_oxygens:
                    glycosidic_bonds += 1

    # For a tetrasaccharide, there should be at least 3 glycosidic bonds connecting the 4 units
    if glycosidic_bonds < 3:
        return False, f"Found {glycosidic_bonds} glycosidic bonds, need at least 3"

    return True, "Contains four monosaccharide units linked via glycosidic bonds"