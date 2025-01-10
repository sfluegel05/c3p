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
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if not nitro_matches:
        return False, "No nitro group found"

    # Ensure that all oxygen atoms are part of nitro groups
    oxygen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    nitro_oxygen_atoms = set()
    for match in nitro_matches:
        nitro_oxygen_atoms.add(match[1])  # Index of oxygen in the nitro group pattern
        nitro_oxygen_atoms.add(match[2])  # Index of the other oxygen in the nitro group pattern

    if len(oxygen_atoms) != len(nitro_oxygen_atoms):
        return False, "Contains oxygen atoms not part of nitro groups"

    # Ensure that all nitrogen atoms are part of nitro groups
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    nitro_nitrogen_atoms = set(match[0] for match in nitro_matches)  # Index of nitrogen in the nitro group pattern

    if len(nitrogen_atoms) != len(nitro_nitrogen_atoms):
        return False, "Contains nitrogen atoms not part of nitro groups"

    return True, "Hydrocarbon with at least one nitro group"