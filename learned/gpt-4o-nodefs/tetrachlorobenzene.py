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

    # Define a set of SMARTS patterns for a benzene ring with at least four attached chlorine atoms, 
    # allowing for connectivity to other rings or groups.
    tetra_pattern_1 = Chem.MolFromSmarts("c1(c(Cl)c(Cl)c(Cl)c(Cl)c1)")
    tetra_pattern_2 = Chem.MolFromSmarts("c1c(Cl)c(Cl)c(Cl)c(Cl)c1")
    biphenyl_pattern = Chem.MolFromSmarts("c1(c(c(Cl)c(Cl)c1)C2=CC=C(Cl)C=C2)")
    
    # Check if the molecule contains any of the patterns for tetrachlorobenzene or derivatives
    if mol.HasSubstructMatch(tetra_pattern_1):
        return True, "Contains tetrachlorobenzene core structure or derivative"
    if mol.HasSubstructMatch(tetra_pattern_2):
        return True, "Contains tetrachlorobenzene core structure or derivative"
    if mol.HasSubstructMatch(biphenyl_pattern):
        return True, "Contains tetrachlorobenzene core structure or derivative"
    
    return False, "Does not contain tetrachlorobenzene core structure or derivative"