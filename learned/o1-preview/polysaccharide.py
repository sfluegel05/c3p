"""
Classifies: CHEBI:18154 polysaccharide
"""
"""
Classifies: polysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdqueries

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

    # Kekulize the molecule (necessary for some substructure searches)
    try:
        Chem.Kekulize(mol)
    except:
        pass  # Kekulization may fail for some molecules

    # Define a pattern for monosaccharide units (pyranose and furanose rings)
    # This is a simplified pattern for a 5 or 6 membered ring with one oxygen and multiple hydroxyl groups
    monosaccharide_pattern = Chem.MolFromSmarts("""
        [C,c]1[C,c][C,c][C,c][C,c][O]1  # 6-membered ring with oxygen (pyranose)
        |
        [C,c]1[C,c][C,c][C,c][O]1       # 5-membered ring with oxygen (furanose)
    """)
    if monosaccharide_pattern is None:
        return False, "Failed to create monosaccharide SMARTS pattern"

    # Find all monosaccharide units
    mono_matches = mol.GetSubstructMatches(monosaccharide_pattern)
    num_monosaccharides = len(mono_matches)

    if num_monosaccharides == 0:
        return False, "No monosaccharide units found"

    if num_monosaccharides <= 10:
        return False, f"Only {num_monosaccharides} monosaccharide units found, need more than 10"

    # Check for glycosidic linkages between monosaccharide units
    # Glycosidic bond pattern: an ether linkage connecting two carbons from different rings
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C]-O-[C]")
    if glycosidic_bond_pattern is None:
        return False, "Failed to create glycosidic bond SMARTS pattern"

    glycosidic_bonds = mol.GetSubstructMatches(glycosidic_bond_pattern)
    num_glycosidic_bonds = len(glycosidic_bonds)

    if num_glycosidic_bonds < (num_monosaccharides - 1):
        return False, f"Insufficient glycosidic linkages: found {num_glycosidic_bonds}, expected at least {num_monosaccharides - 1}"

    # Additional check: Ensure that glycosidic bonds connect different monosaccharide units
    # Build a graph of monosaccharide units connected via glycosidic bonds
    from collections import defaultdict
    mono_atom_indices = [atom_idx for match in mono_matches for atom_idx in match]
    mono_atoms = set(mono_atom_indices)

    glyco_bonds_connecting_monosaccharides = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
            continue
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if bond.GetSmarts() != "O":
            continue
        if (atom1.GetIdx() in mono_atoms) and (atom2.GetIdx() in mono_atoms):
            glyco_bonds_connecting_monosaccharides += 1

    if glyco_bonds_connecting_monosaccharides < (num_monosaccharides - 1):
        return False, f"Insufficient glycosidic bonds connecting monosaccharide units: found {glyco_bonds_connecting_monosaccharides}"

    return True, f"Contains {num_monosaccharides} monosaccharide units connected via glycosidic bonds"