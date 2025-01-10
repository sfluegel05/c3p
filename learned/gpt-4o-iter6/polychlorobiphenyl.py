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
    
    # Define the biphenyl SMARTS pattern
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")
    
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return False, "No biphenyl core structure found"

    # Count chlorine atoms
    chlorine_atoms = [
        atom for atom in mol.GetAtoms()
        if atom.GetSymbol() == "Cl"
    ]
    chlorine_count = len(chlorine_atoms)

    # Ensure biphenyl rings are substituted with chlorine between 2 to 10 times
    if chlorine_count < 2 or chlorine_count > 10:
        return False, f"Found {chlorine_count} chlorine atoms, should be between 2 and 10"
    
    # Using more comprehensive matching for other potential functional groups
    disruptive_groups_smarts = "[OH,N(=O)=O,!Cl]"
    disruptive_groups = Chem.MolFromSmarts(disruptive_groups_smarts)
    if mol.HasSubstructMatch(disruptive_groups):
        return False, "Contains disruptive functional groups like hydroxyls, nitro, or others"

    return True, f"Contains biphenyl structure with {chlorine_count} chlorine atoms attached"