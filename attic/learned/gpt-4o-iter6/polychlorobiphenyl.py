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
    
    # Match biphenyl core structure pattern
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")
    biphenyl_match = mol.GetSubstructMatch(biphenyl_pattern)
    if not biphenyl_match:
        return False, "No biphenyl core structure found"

    # Count chlorine atoms directly attached to the biphenyl core carbons
    chlorine_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "Cl":
            for bonded_atom in atom.GetNeighbors():
                if bonded_atom.GetIdx() in biphenyl_match and bonded_atom.GetSymbol() == "C":
                    chlorine_count += 1
                    break

    # Ensure biphenyl rings are substituted with chlorine between 2 to 10 times
    if chlorine_count < 2 or chlorine_count > 10:
        return False, f"Found {chlorine_count} chlorine atoms on biphenyl rings, should be between 2 and 10"

    return True, f"Contains biphenyl structure with {chlorine_count} chlorine atoms attached"