"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene or its derivative based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrachlorobenzene or its derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES into RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for a benzene ring with four attached chlorine atoms
    tetrachloro_benzene_pattern = Chem.MolFromSmarts("c1cc(Cl)cc(Cl)c1(Cl)Cl")
    
    # Check if the molecule contains a tetrachlorobenzene pattern
    if mol.HasSubstructMatch(tetrachloro_benzene_pattern):
        return True, "Contains tetrachlorobenzene core structure"
    
    return False, "Does not contain tetrachlorobenzene core structure"