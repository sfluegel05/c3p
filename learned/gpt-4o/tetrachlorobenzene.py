"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is defined as a benzene ring carrying four chlorine groups at unspecified positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrachlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for a benzene ring with exactly four chlorines
    tetrachloro_pattern = Chem.MolFromSmarts("c1c(Cl)c(Cl)c(Cl)c(Cl)c1")
    
    # Check if the pattern is present
    if mol.HasSubstructMatch(tetrachloro_pattern):
        return True, "Contains a benzene ring with four chlorine atoms"
    
    return False, "Does not match the tetrachlorobenzene pattern"

# Example Usage:   
# print(is_tetrachlorobenzene("Clc1cc(Cl)c(Cl)cc1Cl"))
# This should return True, "Contains a benzene ring with four chlorine atoms"