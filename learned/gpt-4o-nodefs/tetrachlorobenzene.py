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

    # Define a more comprehensive set of SMARTS patterns for a benzene ring with at least four chlorine atoms.
    # Add more patterns as we identify other structural motifs.
    tetra_chloro_phenol_pattern = Chem.MolFromSmarts("c1(c(O)c(Cl)c(Cl)c(Cl)c(Cl)c1)")
    pentachloro_pattern = Chem.MolFromSmarts("c1(c(Cl)c(Cl)c(Cl)c(Cl)c(Cl)c1)")
    tetrachlorobenzene_pattern = Chem.MolFromSmarts("c1c(Cl)c(Cl)c(Cl)c(Cl)c1")

    # Check if the molecule contains any of the updated patterns for tetrachlorobenzene or derivatives
    if mol.HasSubstructMatch(tetra_chloro_phenol_pattern):
        return True, "Contains tetrachloro phenol structure or derivative"
    if mol.HasSubstructMatch(pentachloro_pattern):
        return True, "Contains pentachloro benzene structure"
    if mol.HasSubstructMatch(tetrachlorobenzene_pattern):
        return True, "Contains tetrachlorobenzene structure"

    return False, "Does not contain tetrachlorobenzene core structure or derivative"