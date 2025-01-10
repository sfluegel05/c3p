"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: CHEBI:82795 tetrachlorobenzene
"""
from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule contains a tetrachlorobenzene moiety based on its SMILES string.
    A tetrachlorobenzene moiety is a benzene ring with exactly four chlorine atoms attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a tetrachlorobenzene moiety, False otherwise
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

    # Count the total number of chlorine atoms in the molecule
    chlorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17)
    if chlorine_count < 4:
        return False, f"Found {chlorine_count} chlorine atoms, need at least 4"

    # Check if the benzene ring has exactly four chlorine atoms attached
    tetrachlorobenzene_pattern = Chem.MolFromSmarts("c1(Cl)c(Cl)c(Cl)c(Cl)cc1")
    if mol.HasSubstructMatch(tetrachlorobenzene_pattern):
        return True, "Contains a benzene ring with exactly four chlorine atoms attached"

    return False, "No benzene ring with exactly four chlorine atoms found"