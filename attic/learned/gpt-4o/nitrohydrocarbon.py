"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon in which one or more of the hydrogens has been replaced
    by nitro groups.

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
    if not nitro_matches:
        return False, "No nitro groups found"

    # Ensure the molecule is predominantly a hydrocarbon
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    non_carbon_non_nitro_count = sum(1 for atom in mol.GetAtoms() 
                                     if atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 7 and atom.GetAtomicNum() != 8)

    # There should be more carbons than other non-nitro/non-oxy atoms
    if carbon_count <= non_carbon_non_nitro_count:
        return False, "Molecule is not predominantly a hydrocarbon"

    return True, "Contains nitro group(s) attached to hydrocarbon backbone"