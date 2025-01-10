"""
Classifies: CHEBI:27096 trichlorobenzene
"""
from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene based on its SMILES string.
    A trichlorobenzene is specifically a benzene ring with three chloro substituents.

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

    # Define a SMARTS pattern for a benzene ring with specifically three Cl atoms
    trichlorobenzene_pattern = Chem.MolFromSmarts("c1c(Cl)cc(Cl)cc1Cl")
    
    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(trichlorobenzene_pattern):
        return True, "Contains a benzene ring with exactly three chloro substituents"
    
    return False, "Does not contain a benzene ring with exactly three chloro substituents"