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
    and without glycosidic connections to other saccharide units.

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

    if not has_aldehyde and not has_ketone:
        return False, "Does not contain an aldehyde or ketone group"

    # Check for multiple hydroxyl groups (-OH) attached to carbons
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4;H2,H1](O)")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, f"Contains {len(hydroxyl_matches)} hydroxyl groups attached to carbons, which is less than 2"

    # Ensure molecule is a single saccharide unit (no glycosidic bonds)
    # Glycosidic bonds are C-O-C links between saccharide units
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C;!R;!$(C(=O))]-O-[C;!R;!$(C(=O))]")
    if mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "Contains glycosidic bonds indicating connection to other saccharide units"

    # Optional: Check for deoxy sugars (allowing for missing hydroxyls)
    # Deoxy sugars are allowed as per definition, so we do not enforce all carbons have hydroxyls

    # Monosaccharides can be cyclic (furanose or pyranose rings) or open-chain forms
    # We accept both forms, so no need to enforce ring structures

    return True, "Molecule satisfies criteria for a monosaccharide"