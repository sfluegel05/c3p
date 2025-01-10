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

    # Look for biphenyl core structure - two phenyl rings directly connected
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return False, "No biphenyl core structure found"

    # Count the number of chlorine atoms attached to the entire molecule
    chlorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "Cl")

    # Ensure biphenyl rings are substituted with chlorine between 2 to 10 times
    if chlorine_count < 2 or chlorine_count > 10:
        return False, f"Found {chlorine_count} chlorine atoms, should be between 2 and 10"

    # Now, verifying that these cholorines are attached to the biphenyl rings
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "Cl":
            # Ensure that this chlorine atom is directly attached to a carbon
            connected_to_biphenyl = any(neighbor.GetSymbol() == "C" for neighbor in atom.GetNeighbors())
            if not connected_to_biphenyl:
                return False, "Chlorine atoms must be attached to biphenyl rings"

    return True, f"Contains biphenyl structure with {chlorine_count} chlorine atoms attached"