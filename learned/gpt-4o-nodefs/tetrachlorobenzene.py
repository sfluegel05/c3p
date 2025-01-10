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

    # Define a more generalized SMARTS pattern for a benzene ring with at least four attached chlorine atoms
    generalized_pattern = Chem.MolFromSmarts("c1c(Cl)c(Cl)c(Cl)c(Cl)c1")
    
    # Check if the molecule contains any general tetrachlorobenzene or derivatives
    if mol.HasSubstructMatch(generalized_pattern):
        return True, "Contains tetrachlorobenzene core structure or derivative"

    return False, "Does not contain tetrachlorobenzene core structure or derivative"