"""
Classifies: CHEBI:35746 fatty aldehyde
"""
from rdkit import Chem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is characterized by a terminal aldehyde group at the end of a long aliphatic chain,
    possibly including unsaturation (double bonds).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to get the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern to identify a terminal aldehyde group
    # Identifying a CHO group with the carbon (C=O) having a single hydrogen (H) attached to it
    aldehyde_pattern = Chem.MolFromSmarts('[CX3H1](=O)')
    
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)

    if not aldehyde_matches:
        return False, "No terminal aldehyde group found"

    # Check for a sufficient aliphatic carbon chain attached to the aldehyde
    def is_valid_aliphatic(atom):
        # Carbon atom not in a ring, part of a continuous chain
        return atom.GetAtomicNum() == 6 and not atom.IsInRing()

    for match in aldehyde_matches:
        aldehyde_carbon_idx = match[0]  

        # Start chain length counting from the neighbor of the aldehyde carbon
        chain_length = 0
        visited = set()

        def traverse_chain(atom_idx, prev_idx=None):
            nonlocal chain_length
            if atom_idx not in visited:
                visited.add(atom_idx)
                atom = mol.GetAtomWithIdx(atom_idx)
                if is_valid_aliphatic(atom):
                    chain_length += 1
                    # Traverse neighbors excluding the previous atom
                    for neighbor in atom.GetNeighbors():
                        next_idx = neighbor.GetIdx()
                        if next_idx != prev_idx:
                            traverse_chain(next_idx, atom_idx)
        
        # Ensure the aldehyde carbon has exactly one sp3 carbon neighbor, indicating terminality
        carbon_neighbors = [nbr for nbr in mol.GetAtomWithIdx(aldehyde_carbon_idx).GetNeighbors() if nbr.GetAtomicNum() == 6]
        
        if len(carbon_neighbors) != 1:
            return False, "Aldehyde is not at the terminal end of a chain"

        # Traverse from that single neighbor
        traverse_chain(carbon_neighbors[0].GetIdx(), aldehyde_carbon_idx)

        # Ensure the chain is long and complex enough to represent a typical fatty aldehyde
        if chain_length >= 6:
            return True, "Valid fatty aldehyde: Terminal aldehyde group with a suitable aliphatic chain"
    
    return False, "Carbon chain too short for typical fatty aldehyde"