"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl based on its SMILES string.
    A polychlorobiphenyl is defined as a biphenyl compound containing between 2 and 10 chlorine atoms
    attached to the two benzene rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polychlorobiphenyl, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify biphenyl structure: two phenyl rings joined via a single bond
    # SMARTS pattern for biphenyl structure:
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return False, "No biphenyl structure found"

    # Count chlorine atoms attached to the biphenyl rings
    chlorine_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17]
    connected_to_biphenyl = sum(
        any(neighbor.HasSubstructMatch(Chem.MolFromSmarts("c1ccccc1")) or 
            neighbor.HasSubstructMatch(Chem.MolFromSmarts("c2ccccc2"))
            for neighbor in chlorine.GetNeighbors()) 
        for chlorine in chlorine_atoms
    )

    # Check the count of chlorines
    if connected_to_biphenyl < 2 or connected_to_biphenyl > 10:
        return False, f"Contains {connected_to_biphenyl} chlorine atoms connected to biphenyl, needs between 2 and 10"

    return True, "Contains biphenyl structure with acceptable number of chlorine atoms (between 2 and 10)"