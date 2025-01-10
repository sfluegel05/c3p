"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Check for two connected benzene rings (biphenyl)
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return False, "No biphenyl structure found"

    # Count chlorine atoms attached to carbons in benzene rings
    chlorine_pattern = Chem.MolFromSmarts("[cH]Cl")
    chlorine_matches = mol.GetSubstructMatches(chlorine_pattern)
    num_chlorines = len(chlorine_matches)

    # Check if the number of chlorines is between 2 and 10
    if num_chlorines < 2 or num_chlorines > 10:
        return False, f"Found {num_chlorines} chlorine atoms, should be between 2 and 10"

    return True, f"Contains biphenyl structure with {num_chlorines} chlorine atoms attached"