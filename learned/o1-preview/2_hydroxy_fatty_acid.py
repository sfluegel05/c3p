"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: CHEBI:15617 2-hydroxy fatty acid
"""
from rdkit import Chem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid is any fatty acid with a hydroxy functional group in the alpha- or 2-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of only C, H, O atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6, 8]:
            return False, "Molecule contains atoms other than C, H, and O"

    # Check for exactly one carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, "Molecule does not have exactly one carboxylic acid group"
    carboxyl_c_index = carboxylic_acid_matches[0][0]

    # Check for rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring structures, atypical for fatty acids"

    # Check for exactly one hydroxy group (-OH)
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) != 1:
        return False, "Molecule does not have exactly one hydroxy group"

    # Check that the hydroxy group is on the alpha carbon adjacent to carboxyl carbon
    carboxyl_c_atom = mol.GetAtomWithIdx(carboxyl_c_index)
    alpha_c_atom = None
    for neighbor in carboxyl_c_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # Carbon atom
            alpha_c_atom = neighbor
            break
    if alpha_c_atom is None:
        return False, "No alpha carbon found adjacent to carboxyl carbon"
    alpha_c_index = alpha_c_atom.GetIdx()

    # Check that the hydroxy group is attached to the alpha carbon
    hydroxy_on_alpha = False
    for neighbor in alpha_c_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
            # Oxygen atom with single bond (hydroxy group)
            hydroxy_on_alpha = True
            break
    if not hydroxy_on_alpha:
        return False, "Hydroxy group is not on the alpha carbon"

    # Ensure no other functional groups (only C, H, O)
    # We've already checked atom types and the number of hydroxy groups

    # Optionally, check molecule is a fatty acid (e.g., chain length)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, f"Molecule has {c_count} carbon atoms, less than typical fatty acids"

    return True, "Molecule is a 2-hydroxy fatty acid with hydroxy group at alpha position"