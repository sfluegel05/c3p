"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: polyprenol
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol is any prenol possessing the general formula H-[CH2C(Me)=CHCH2]nOH
    in which the carbon skeleton is composed of more than one isoprene units.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for isoprene unit in polyprenol
    # Matches the pattern: C=C-C-C
    isoprene_unit_smarts = "[CH2]=[CH]-[CH2]-[CH2]"
    isoprene_unit = Chem.MolFromSmarts(isoprene_unit_smarts)
    
    if isoprene_unit is None:
        return False, "Failed to create isoprene unit pattern"
    
    # Find all matches of the isoprene unit
    isoprene_matches = mol.GetSubstructMatches(isoprene_unit)
    num_isoprene_units = len(isoprene_matches)
    
    if num_isoprene_units <= 1:
        return False, f"Found {num_isoprene_units} isoprene units, need more than 1"
    
    # Check if the isoprene units are connected head-to-tail
    # This can be complex, so we'll assume that if the units are connected linearly, it's acceptable
    # Get all atoms involved in isoprene units
    isoprene_atoms = set()
    for match in isoprene_matches:
        isoprene_atoms.update(match)
    
    # Check if the isoprene atoms form a continuous chain
    # Get the shortest path between the first and last isoprene atom
    path = Chem.rdmolops.GetShortestPath(mol, min(isoprene_atoms), max(isoprene_atoms))
    if not set(path).issubset(isoprene_atoms):
        return False, "Isoprene units are not connected linearly"
    
    # Check for terminal hydroxyl group (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No terminal hydroxyl group found"
    
    # Verify that the hydroxyl group is at the end of the isoprene chain
    terminal_carbons = [atom_idx for atom_idx in isoprene_atoms if mol.GetAtomWithIdx(atom_idx).GetDegree() == 1]
    hydroxyl_connected = False
    for match in hydroxyl_matches:
        hydroxyl_carbon_idx = match[0]
        if hydroxyl_carbon_idx in terminal_carbons:
            hydroxyl_connected = True
            break
    if not hydroxyl_connected:
        return False, "Hydroxyl group is not connected to the end of the isoprene chain"
    
    return True, f"Contains {num_isoprene_units} isoprene units connected head-to-tail with terminal hydroxyl group"