"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
from rdkit import Chem

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl based on its SMILES string.
    A polychlorobiphenyl is a biphenyl compound containing between 2 and 10 chlorine atoms
    attached to the two benzene rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polychlorobiphenyl, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved pattern for biphenyl structure with chlorine atoms
    pattern = Chem.MolFromSmarts("c1c(Cl)cccc1-c2c(Cl)cccc2")
    if not mol.HasSubstructMatch(pattern):
        return False, "No biphenyl structure with appropriate chlorine atoms found"

    # Count total chlorine atoms
    chlorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "Cl")
    
    # Check if the number of chlorines is between 2 and 10
    if chlorine_count < 2 or chlorine_count > 10:
        return False, f"Found {chlorine_count} chlorine atoms, should be between 2 and 10"

    return True, f"Contains biphenyl structure with {chlorine_count} chlorine atoms attached"