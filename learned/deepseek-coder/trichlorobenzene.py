"""
Classifies: CHEBI:27096 trichlorobenzene
"""
"""
Classifies: CHEBI:36043 trichlorobenzene
"""
from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene based on its SMILES string.
    A trichlorobenzene is a benzene ring with exactly three chlorine substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trichlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for benzene ring pattern
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene_pattern):
        return False, "No benzene ring found"

    # Count the number of chlorine atoms attached to the benzene ring
    chlorine_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 17:  # Chlorine has atomic number 17
            # Check if the chlorine is attached to a benzene ring
            for neighbor in atom.GetNeighbors():
                if neighbor.IsInRing() and neighbor.GetAtomicNum() == 6:  # Carbon in a ring
                    chlorine_count += 1
                    break

    if chlorine_count == 3:
        return True, "Contains a benzene ring with exactly three chlorine substituents"
    else:
        return False, f"Found {chlorine_count} chlorine substituents on benzene ring, need exactly 3"