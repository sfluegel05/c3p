"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
"""
Classifies: omega-hydroxy fatty acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from collections import deque

def get_longest_carbon_chain(mol, start_idx):
    """Find the longest continuous carbon chain from a starting point"""
    visited = set()
    longest_path = []
    
    def dfs(atom_idx, current_path):
        nonlocal longest_path
        visited.add(atom_idx)
        current_path.append(atom_idx)
        
        if len(current_path) > len(longest_path):
            longest_path = current_path.copy()
            
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                dfs(neighbor.GetIdx(), current_path)
        current_path.pop()
        visited.remove(atom_idx)
    
    dfs(start_idx, [])
    return longest_path

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
    
    # Look for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[OX2H][CX3](=[OX1])")
    cooh_matches = mol.GetSubstructMatches(carboxylic_pattern)
    
    if not cooh_matches:
        return False, "No carboxylic acid group found"
    if len(cooh_matches) > 1:
        return False, "Multiple carboxylic acid groups found"
    
    # Get the carboxylic acid carbon
    cooh_carbon_idx = cooh_matches[0][1]
    
    # Find terminal hydroxyl groups
    terminal_oh_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    oh_matches = mol.GetSubstructMatches(terminal_oh_pattern)
    
    if not oh_matches:
        return False, "No hydroxyl groups found"
    
    # Find valid terminal hydroxyl groups
    valid_terminal_oh = []
    for match in oh_matches:
        carbon_idx = match[0]
        carbon = mol.GetAtomWithIdx(carbon_idx)
        
        # Check if carbon is terminal (only one carbon neighbor)
        carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            valid_terminal_oh.append(carbon_idx)
    
    if not valid_terminal_oh:
        return False, "No terminal hydroxyl groups found"
    
    # Check ring characteristics
    ring_info = mol.GetRingInfo()
    ring_atoms = set()
    for ring in ring_info.AtomRings():
        ring_atoms.update(ring)
        if len(ring) > 3:  # Allow only small (epoxy) rings
            return False, "Contains large ring structure"
    
    # Find the longest carbon chain from COOH
    main_chain = get_longest_carbon_chain(mol, cooh_carbon_idx)
    if len(main_chain) < 3:
        return False, "Carbon chain too short"
    
    # Check if any terminal OH is at the end of the longest chain
    valid_omega_oh = False
    for oh_carbon in valid_terminal_oh:
        if oh_carbon in main_chain and main_chain.index(oh_carbon) in [0, len(main_chain)-1]:
            valid_omega_oh = True
            break
    
    if not valid_omega_oh:
        return False, "No omega hydroxyl group found"
    
    # Additional structural checks
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            # Check for excessive branching
            carbon_neighbors = len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6])
            if carbon_neighbors > 2 and atom.GetIdx() not in ring_atoms:
                return False, "Structure too branched"
        elif atom.GetAtomicNum() not in [1, 6, 8]:
            return False, "Contains elements other than C, H, O"
    
    return True, "Contains terminal (omega) hydroxyl and carboxylic acid in a linear structure"