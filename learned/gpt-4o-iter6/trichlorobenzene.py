"""
Classifies: CHEBI:27096 trichlorobenzene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene based on its SMILES string.
    A trichlorobenzene is a member of the chlorobenzenes class with three chloro substituents on the benzene ring.

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

    # SMARTS pattern for a benzene ring with three chlorine substituents
    trichlorobenzene_pattern = Chem.MolFromSmarts("c1c(Cl)c(Cl)c(Cl)ccc1")
    if mol.HasSubstructMatch(trichlorobenzene_pattern):
        return True, "Contains a benzene ring with three chloro substituents"

    return False, "Benzene ring with three chloro substituents not found"