"""
Classifies: CHEBI:27096 trichlorobenzene
"""
from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene based on its SMILES string.
    A trichlorobenzene has a benzene ring with exactly three chlorine atoms attached to it.

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

    # Check for benzene ring
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene_pattern):
        return False, "No benzene ring found"
    
    # Get all chlorine atoms
    chlorine_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17]
    
    # Check if three chlorines are directly attached to the benzene ring
    benzene_match = mol.GetSubstructMatch(benzene_pattern)
    chlorine_count = 0
    for atom in chlorine_atoms:
        for neighbor in atom.GetNeighbors():
             if neighbor.GetIdx() in benzene_match:
                chlorine_count += 1
                
    if chlorine_count != 3:
            return False, f"Found {chlorine_count} chlorine atoms attached to benzene, need exactly 3"
        
    return True, "Contains a benzene ring with exactly 3 chlorine atoms attached"