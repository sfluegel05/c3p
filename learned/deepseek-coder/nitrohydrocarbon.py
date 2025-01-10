"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
"""
Classifies: CHEBI:51050 nitrohydrocarbon
"""
from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon in which one or more hydrogens have been replaced by nitro groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for elements other than C, H, N, O (nitro group is allowed)
    allowed_elements = {6, 1, 7, 8}  # C, H, N, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_elements:
            return False, f"Contains element {atom.GetSymbol()}, which is not allowed in nitrohydrocarbons"

    # Check for at least one nitro group
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    if not mol.HasSubstructMatch(nitro_pattern):
        return False, "No nitro group found"

    # Ensure that the molecule is a hydrocarbon (only C and H, except for nitro groups)
    # Count the number of non-C, non-H atoms (excluding nitro groups)
    non_CH_atoms = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in {6, 1}:  # Not C or H
            non_CH_atoms += 1

    # Subtract the number of nitro groups (each nitro group has 1 N and 2 O atoms)
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    nitro_group_count = len(nitro_matches)
    non_CH_atoms -= 3 * nitro_group_count  # Each nitro group contributes 3 non-CH atoms (1 N + 2 O)

    if non_CH_atoms > 0:
        return False, "Contains non-hydrocarbon atoms other than nitro groups"

    return True, "Hydrocarbon with at least one nitro group"