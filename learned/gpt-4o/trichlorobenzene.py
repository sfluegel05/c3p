"""
Classifies: CHEBI:27096 trichlorobenzene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene based on its SMILES string.
    A trichlorobenzene is a benzene ring with three chlorine substituents.

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

    # Define benzene pattern
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene_pattern):
        return False, "No benzene ring found"
    
    # Check for chlorine atoms directly attached to the benzene ring
    chlorine_pattern = Chem.MolFromSmarts("c[Cl]")
    chlorine_matches = mol.GetSubstructMatches(chlorine_pattern)

    # Count chlorine atoms attached to benzene
    chlorine_count = len(chlorine_matches)
    
    # Check if there are exactly three chlorine atoms
    if chlorine_count != 3:
        return False, f"Found {chlorine_count} chlorine substituents, needs exactly 3 attached to benzene"

    return True, "Contains benzene ring with three chlorine substituents"