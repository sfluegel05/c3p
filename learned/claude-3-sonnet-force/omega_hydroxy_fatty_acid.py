"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
"""
Classifies: omega-hydroxy fatty acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from collections import deque

def find_path_to_cooh(mol, start_atom_idx, cooh_carbon_idx):
    """Helper function to find shortest path between two atoms"""
    visited = set()
    queue = deque([(start_atom_idx, [start_atom_idx])])
    
    while queue:
        current_idx, path = queue.popleft()
        if current_idx == cooh_carbon_idx:
            return path
            
        atom = mol.GetAtomWithIdx(current_idx)
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx not in visited:
                visited.add(n_idx)
                new_path = list(path)
                new_path.append(n_idx)
                queue.append((n_idx, new_path))
    return None

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid groups
    carboxylic_pattern = Chem.MolFromSmarts("[OX2H][CX3](=[OX1])")
    cooh_matches = mol.GetSubstructMatches(carboxylic_pattern)
    
    if not cooh_matches:
        return False, "No carboxylic acid group found"
    if len(cooh_matches) > 1:
        return False, "Multiple carboxylic acid groups found"
    
    # Look for amide bonds (to exclude peptides)
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    if mol.HasSubstructMatch(amide_pattern):
        return False, "Contains amide bonds (peptide-like structure)"
    
    # Get the carboxylic acid carbon
    cooh_carbon_idx = cooh_matches[0][1]
    
    # Find hydroxyl groups (excluding the one in COOH)
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    oh_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    if not oh_matches:
        return False, "No suitable hydroxyl groups found"
    
    # Check carbon chain characteristics
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 3:
        return False, "Carbon chain too short"
    
    # Find the terminal carbon with hydroxyl that's furthest from COOH
    max_path_length = 0
    valid_terminal_oh = False
    
    for match in oh_matches:
        carbon_idx = match[0]  # Carbon attached to OH
        atom = mol.GetAtomWithIdx(carbon_idx)
        
        # Check if this carbon is terminal (only one carbon neighbor)
        carbon_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            continue
            
        # Find path to COOH
        path = find_path_to_cooh(mol, carbon_idx, cooh_carbon_idx)
        if path and len(path) > max_path_length:
            max_path_length = len(path)
            valid_terminal_oh = True
    
    if not valid_terminal_oh:
        return False, "No valid terminal hydroxyl group found"
    
    # Verify chain characteristics
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 1:  # Allow up to 1 ring (for epoxy groups)
        return False, "Too many ring structures"
    
    # Count branching points
    branching_points = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            carbon_neighbors = len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6])
            if carbon_neighbors > 2:
                branching_points += 1
    
    if branching_points > 2:  # Allow limited branching
        return False, "Structure too branched"
        
    return True, "Contains terminal (omega) hydroxyl and carboxylic acid in a primarily linear structure"