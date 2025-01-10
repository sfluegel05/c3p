"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: polyprenol
"""
from rdkit import Chem

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
    
    # Define SMARTS pattern for the isoprene repeating unit
    # This pattern matches [CH2]-C(=C)-C(-CH3)-[CH2] with any stereochemistry
    isoprene_unit_smarts = "[CH2]-[C](=C)-[C]-[CH2]"
    isoprene_unit = Chem.MolFromSmarts(isoprene_unit_smarts)
    if isoprene_unit is None:
        return False, "Failed to create isoprene unit pattern"
    
    # Find all matches of the isoprene unit
    isoprene_matches = mol.GetSubstructMatches(isoprene_unit, useChirality=False)
    num_isoprene_units = len(isoprene_matches)
    
    if num_isoprene_units <= 1:
        return False, f"Found {num_isoprene_units} isoprene units, need more than 1"
    
    # Check that isoprene units are connected head-to-tail
    # Create a list of atom indices for all isoprene units
    isoprene_atom_indices = [set(match) for match in isoprene_matches]
    
    # Sort units based on their position in the molecule
    isoprene_atom_indices.sort(key=lambda x: min(x))
    
    # Verify connectivity between consecutive isoprene units
    is_connected = True
    for i in range(len(isoprene_atom_indices) - 1):
        unit1 = isoprene_atom_indices[i]
        unit2 = isoprene_atom_indices[i + 1]
        connected = False
        for idx1 in unit1:
            atom1 = mol.GetAtomWithIdx(idx1)
            for neighbor in atom1.GetNeighbors():
                if neighbor.GetIdx() in unit2:
                    connected = True
                    break
            if connected:
                break
        if not connected:
            is_connected = False
            break
    if not is_connected:
        return False, "Isoprene units are not connected head-to-tail"
    
    # Check for terminal hydroxyl group (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No terminal hydroxyl group found"
    
    # Verify that the hydroxyl group is at the end of the isoprene chain
    # Get the indices of all atoms in the isoprene units
    isoprene_atoms = set()
    for unit in isoprene_atom_indices:
        isoprene_atoms.update(unit)
    
    # Check if hydroxyl oxygen is connected to an atom in the isoprene chain
    hydroxyl_connected = False
    for match in hydroxyl_matches:
        hydroxyl_oxygen_idx = match[0]
        hydroxyl_oxygen = mol.GetAtomWithIdx(hydroxyl_oxygen_idx)
        for neighbor in hydroxyl_oxygen.GetNeighbors():
            if neighbor.GetIdx() in isoprene_atoms:
                hydroxyl_connected = True
                break
        if hydroxyl_connected:
            break
    if not hydroxyl_connected:
        return False, "Hydroxyl group is not connected to the isoprene chain"
    
    return True, f"Contains {num_isoprene_units} isoprene units connected head-to-tail with terminal hydroxyl group"