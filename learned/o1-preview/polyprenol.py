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
    A polyprenol is a prenol with more than one isoprene units connected head-to-tail, ending with an -OH group.
    
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
    
    # Remove hydrogens for simplicity
    mol = Chem.RemoveHs(mol)
    
    # Define SMARTS pattern for isoprene unit connected head-to-tail
    # Isoprene unit: CH2=C(CH3)-CH=CH2
    isoprene_pattern = Chem.MolFromSmarts("C(=C(C))C=C")
    if isoprene_pattern is None:
        return False, "Failed to create isoprene unit pattern"
    
    # Find all isoprene units
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    num_isoprene_units = len(isoprene_matches)
    
    if num_isoprene_units <= 1:
        return False, f"Found {num_isoprene_units} isoprene units, need more than 1"
    
    # Verify isoprene units are connected head-to-tail
    # Create a list of isoprene unit indices
    isoprene_units = []
    for match in isoprene_matches:
        isoprene_units.append(match)
    
    # Sort the isoprene units based on atom indices to reflect their order in the molecule
    isoprene_units.sort(key=lambda x: min(x))
    
    # Check connectivity between isoprene units
    is_connected = True
    for i in range(len(isoprene_units) - 1):
        current_unit = isoprene_units[i]
        next_unit = isoprene_units[i + 1]
        connected = False
        for atom_idx in current_unit:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in next_unit:
                    connected = True
                    break
            if connected:
                break
        if not connected:
            is_connected = False
            break
    
    if not is_connected:
        return False, "Isoprene units are not connected head-to-tail"
    
    # Check for terminal hydroxyl group (-CH2OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4H2][OX2H]")  # Primary alcohol
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    if not hydroxyl_matches:
        return False, "No terminal -OH group found"
    
    # Verify the hydroxyl group is at the end of the chain
    # Find the carbon chain length
    chains = mol.GetSubstructMatches(Chem.MolFromSmarts("C[*]"))
    chain_atom_indices = [atom_idx for match in chains for atom_idx in match]
    max_chain_length = len(set(chain_atom_indices))
    
    # Check if the hydroxyl group is connected to one end of the longest carbon chain
    terminal_hydroxyl = False
    for match in hydroxyl_matches:
        hydroxyl_carbon_idx = match[0]
        path_lengths = Chem.rdmolops.GetShortestPath(mol, hydroxyl_carbon_idx, isoprene_units[0][0])
        if len(path_lengths) == max_chain_length - 1:
            terminal_hydroxyl = True
            break
        path_lengths = Chem.rdmolops.GetShortestPath(mol, hydroxyl_carbon_idx, isoprene_units[-1][-1])
        if len(path_lengths) == max_chain_length - 1:
            terminal_hydroxyl = True
            break
    
    if not terminal_hydroxyl:
        return False, "Hydroxyl group is not at the end of the carbon chain"
    
    return True, f"Contains {num_isoprene_units} isoprene units connected head-to-tail with terminal -OH group"