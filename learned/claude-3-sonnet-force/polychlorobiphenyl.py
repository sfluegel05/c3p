"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
"""
Classifies: CHEBI:35448 polychlorobiphenyl
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl based on its SMILES string.
    A polychlorobiphenyl is a biphenyl compound containing between 2 and 10 chlorine atoms attached to the two benzene rings.

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

    # Check for biphenyl core
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c1ccccc1")
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return False, "No biphenyl core found"

    # Count number of chlorine atoms
    num_chlorines = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17)
    if num_chlorines < 2 or num_chlorines > 10:
        return False, f"Found {num_chlorines} chlorine atoms, should be between 2 and 10"

    # Check that chlorines are attached to the benzene rings
    chlorine_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17]
    for cl_atom in chlorine_atoms:
        if not any(neighbor_atom.IsInRingSize(6) for neighbor_atom in cl_atom.GetNeighbors()):
            return False, "Chlorine atom not attached to a benzene ring"

    return True, "Biphenyl compound with 2-10 chlorine atoms attached to the benzene rings"