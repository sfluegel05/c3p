"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: fatty aldehyde
"""
from rdkit import Chem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is an aldehyde arising from reduction of the carboxylic acid group of a fatty acid,
    having a carbonyl group at one end of the carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify aldehyde groups (-CHO)
    # Aldehyde carbon: sp2 carbon (not in ring) double-bonded to oxygen, single-bonded to hydrogen
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not aldehyde_matches:
        return False, "No aldehyde group found"
    
    # Check each aldehyde group
    for match in aldehyde_matches:
        aldehyde_carbon_idx = match[0]
        aldehyde_carbon = mol.GetAtomWithIdx(aldehyde_carbon_idx)
        
        # Get neighbors of aldehyde carbon
        neighbors = aldehyde_carbon.GetNeighbors()
        
        # Count number of carbon neighbors (excluding the carbonyl oxygen)
        carbon_neighbors = [nbr for nbr in neighbors if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            continue  # Not a terminal aldehyde group
        
        # Ensure aldehyde carbon is not in a ring
        if aldehyde_carbon.IsInRing():
            continue  # Aldehyde carbon is in a ring, skip
        
        # Start from the carbon neighbor, traverse the chain
        chain_atom = carbon_neighbors[0]
        visited = set()
        atoms_queue = [chain_atom]
        is_aliphatic_chain = True
        while atoms_queue:
            atom = atoms_queue.pop(0)
            atom_idx = atom.GetIdx()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            
            # Check for aromaticity and rings
            if atom.GetIsAromatic() or atom.IsInRing():
                is_aliphatic_chain = False
                break
            
            # Allow carbon and certain heteroatoms (e.g., oxygen in hydroxyl groups)
            atomic_num = atom.GetAtomicNum()
            if atomic_num not in (6, 8, 1):  # Carbon, Oxygen, Hydrogen
                is_aliphatic_chain = False
                break
            
            # Add neighbors to queue (excluding the aldehyde carbon)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx not in visited and nbr_idx != aldehyde_carbon_idx:
                    atoms_queue.append(nbr)
        if not is_aliphatic_chain:
            continue  # Not an aliphatic chain, check next aldehyde
        else:
            return True, "Contains terminal aldehyde group attached to an aliphatic chain"
    return False, "No terminal aldehyde group attached to an aliphatic chain found"