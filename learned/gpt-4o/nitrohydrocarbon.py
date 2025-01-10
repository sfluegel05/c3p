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
    
    # Define SMARTS pattern for generic nitro group (positive nitrogen bound to two oxygens)
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    
    # Search for nitro groups attached to any carbon atom
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            neighbors = atom.GetNeighbors()
            for neighbor in neighbors:
                if neighbor.GetSmarts() == "[N+](=O)[O-]":
                    return True, "Contains nitro group(s) attached to hydrocarbon"

    # Ensure the molecule is predominantly a hydrocarbon
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    hydrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    
    if carbon_count > 0 and hydrogen_count > 0:
        return False, "No nitro groups attached to carbon found but is a hydrocarbon"
    
    return False, "Molecule is not primarily a hydrocarbon"