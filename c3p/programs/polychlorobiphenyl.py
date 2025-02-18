"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
"""
Classifies: polychlorobiphenyl
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl based on its SMILES string.
    A polychlorobiphenyl has 2-10 chlorine atoms attached to a biphenyl core.

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

    # Look for biphenyl core (two benzene rings connected)
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c1ccccc1")
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return False, "No biphenyl core found"

    # Count chlorine atoms
    cl_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17]
    n_chlorines = len(cl_atoms)
    
    if n_chlorines < 2:
        return False, f"Too few chlorines ({n_chlorines}), minimum is 2"
    if n_chlorines > 10:
        return False, f"Too many chlorines ({n_chlorines}), maximum is 10"
    
    # Verify chlorines are attached to aromatic carbons
    for cl in cl_atoms:
        neighbors = cl.GetNeighbors()
        if len(neighbors) != 1:  # Chlorine should have exactly one neighbor
            return False, "Invalid chlorine bonding"
        neighbor = neighbors[0]
        if not (neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic()):
            return False, "Chlorine not attached to aromatic carbon"
    
    # Get the two benzene rings from biphenyl core
    matches = mol.GetSubstructMatches(biphenyl_pattern)
    if len(matches) > 0:
        ring_atoms = set(matches[0])  # atoms in the biphenyl core
        
        # Check if all chlorines are attached to the rings
        for cl in cl_atoms:
            neighbor = cl.GetNeighbors()[0]
            if neighbor.GetIdx() not in ring_atoms:
                return False, "Chlorine attached outside biphenyl core"
    
    return True, f"Valid polychlorobiphenyl with {n_chlorines} chlorine atoms"