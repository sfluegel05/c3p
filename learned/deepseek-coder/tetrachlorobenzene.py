"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: CHEBI:82795 tetrachlorobenzene
"""
from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is a benzene ring with exactly four chlorine atoms attached.

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

    # Look for benzene ring pattern
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene_pattern):
        return False, "No benzene ring found"

    # Count the total number of chlorine atoms in the molecule
    chlorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17)
    if chlorine_count != 4:
        return False, f"Found {chlorine_count} chlorine atoms, need exactly 4"

    # Ensure the molecule is primarily a benzene ring with four chlorine atoms
    # by checking that the molecule has exactly 10 atoms (6 carbons + 4 chlorines)
    if mol.GetNumAtoms() != 10:
        return False, "Molecule is not a simple tetrachlorobenzene"

    return True, "Contains a benzene ring with exactly four chlorine atoms attached"