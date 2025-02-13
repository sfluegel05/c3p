"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: CHEBI:17522 alditol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is an acyclic polyol with the general formula HOCH2[CH(OH)]nCH2OH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is acyclic (no rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, but alditol must be acyclic"

    # Count the number of hydroxyl groups (OH)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
    
    # Count the number of carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Check if the molecule has the general formula HOCH2[CH(OH)]nCH2OH
    # This means:
    # - The molecule should have at least 2 hydroxyl groups on the terminal carbons
    # - Each internal carbon should have at least one hydroxyl group
    # - The molecule should have at least 3 carbons (n >= 1)

    # Check terminal carbons (should have 2 hydroxyl groups each)
    terminal_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetDegree() == 1]
    if len(terminal_carbons) != 2:
        return False, "Molecule does not have exactly two terminal carbons"

    for terminal_carbon in terminal_carbons:
        neighbors = terminal_carbon.GetNeighbors()
        hydroxyl_neighbors = [neighbor for neighbor in neighbors if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1]
        if len(hydroxyl_neighbors) != 2:
            return False, "Terminal carbons do not have exactly two hydroxyl groups"

    # Check internal carbons (should have at least one hydroxyl group)
    internal_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetDegree() > 1]
    for internal_carbon in internal_carbons:
        neighbors = internal_carbon.GetNeighbors()
        hydroxyl_neighbors = [neighbor for neighbor in neighbors if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1]
        if len(hydroxyl_neighbors) < 1:
            return False, "Internal carbons do not have at least one hydroxyl group"

    # Check the general formula HOCH2[CH(OH)]nCH2OH
    # The number of hydroxyl groups should be equal to the number of carbons + 1
    if hydroxyl_count != carbon_count + 1:
        return False, f"Hydroxyl count ({hydroxyl_count}) does not match the general formula HOCH2[CH(OH)]nCH2OH"

    return True, "Molecule is an acyclic polyol with the general formula HOCH2[CH(OH)]nCH2OH"