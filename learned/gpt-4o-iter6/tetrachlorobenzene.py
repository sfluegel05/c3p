"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is defined as a benzene ring carrying exactly four chloro groups at unspecified positions.

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

    # Define tetrachlorobenzene SMARTS pattern: benzene with 4 chlorines
    tetrachlorobenzene_pattern = Chem.MolFromSmarts("c1(c(Cl)c(Cl)c(Cl)c(Cl)c1)")

    # Search for the pattern in the molecule
    if mol.HasSubstructMatch(tetrachlorobenzene_pattern):
        return True, "Benzene ring with exactly four chloro groups found"

    return False, "No benzene ring with exactly four chloro groups found"