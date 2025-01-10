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

    # Ensure the molecule is primarily a hydrocarbon (composed of carbon and hydrogen)
    # Count the number of carbon and hydrogen atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    total_atoms = mol.GetNumAtoms()

    # Ensure that the majority of atoms are carbon and hydrogen
    if (c_count + h_count) / total_atoms < 0.8:
        return False, "Molecule is not primarily a hydrocarbon"

    # Ensure that the number of carbon atoms is greater than the number of nitro groups
    if c_count <= len(nitro_matches):
        return False, "Molecule is not primarily a hydrocarbon"

    return True, "Hydrocarbon with at least one nitro group"