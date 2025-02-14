"""
Classifies: CHEBI:18133 hexose
"""
"""
Classifies: CHEBI:18133 hexose
"""

from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is any six-carbon monosaccharide which in its linear form contains either an aldehyde group at position 1 (aldohexose) or a ketone group at position 2 (ketohexose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, f"Molecule has {c_count} carbon atoms, expected 6"

    # Check for glycosidic bonds (Oxygen atom connected to two carbons not in a ring)
    glycosidic_bonds = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen atom
            neighbor_carbons = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
            if len(neighbor_carbons) == 2:
                # Oxygen bonded to two carbons
                # Check if oxygen is in a ring
                if not atom.IsInRing():
                    glycosidic_bonds += 1

    if glycosidic_bonds > 0:
        return False, f"Molecule appears to have glycosidic bonds ({glycosidic_bonds} detected), may not be a monosaccharide"

    # Check for aldehyde group (aldohexose) or ketone group (ketohexose)
    # Aldehyde pattern: carbonyl group at the end of carbon chain
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)

    # Ketone pattern: carbonyl group within carbon chain
    ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)

    # In cyclic forms, aldehyde/ketone groups form hemiacetals/hemiketals, so we check for those as well
    hemiacetal_pattern = Chem.MolFromSmarts("[OX2H][CX4H][OX2][CX3]=[OX1]")
    hemiacetal_matches = mol.GetSubstructMatches(hemiacetal_pattern)
    hemiketal_pattern = Chem.MolFromSmarts("[OX2H][CX4H][OX2][CX4H][OX2][CX3]=[OX1]")
    hemiketal_matches = mol.GetSubstructMatches(hemiketal_pattern)

    if not (aldehyde_matches or ketone_matches or hemiacetal_matches or hemiketal_matches):
        return False, "No aldehyde, ketone, hemiacetal, or hemiketal group found"

    # All checks passed
    return True, "Contains 6 carbons, no glycosidic bonds, and contains an aldehyde, ketone, hemiacetal, or hemiketal group"