"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
"""
Classifies: CHEBI:60912 polychlorobiphenyl
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl based on its SMILES string.
    A polychlorobiphenyl is a biphenyl compound containing between 2 and 10 chlorine atoms.

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

    # Check for biphenyl structure (two benzene rings connected by a single bond)
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return False, "No biphenyl structure found"

    # Count the number of chlorine atoms
    chlorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17)
    
    # Check if the number of chlorine atoms is between 2 and 10
    if chlorine_count < 2 or chlorine_count > 10:
        return False, f"Found {chlorine_count} chlorine atoms, need between 2 and 10"

    return True, f"Contains biphenyl structure with {chlorine_count} chlorine atoms"