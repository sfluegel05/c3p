"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon where one or more hydrogens have been replaced by nitro groups.

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

    # Define SMARTS pattern for nitro group
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")

    # Search for nitro groups
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    nitro_count = len(nitro_matches)
    if nitro_count == 0:
        return False, "No nitro groups found"

    # Count atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    hydrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    other_atom_count = sum(1 for atom in mol.GetAtoms()
                           if atom.GetAtomicNum() not in {1, 6, 7, 8})

    # Ensure the molecule has primarily a hydrocarbon character
    if other_atom_count > 0:
        return False, f"Contains {other_atom_count} other types of atoms beyond hydrocarbon and nitro group structures"

    if carbon_count <= 0:
        return False, "Molecule lacks sufficient hydrocarbon backbone"

    # Check if hydrocarbons form the backbone structure
    backbone_ratio = carbon_count / (carbon_count + hydrogen_count)
    if backbone_ratio < 0.5:
        return False, "Hydrophobic carbon backbone too small relative to hydrophilic components"

    return True, "Contains nitro group(s) attached to hydrocarbon backbone as expected for nitrohydrocarbon"