"""
Classifies: CHEBI:35746 fatty aldehyde
"""
from rdkit import Chem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde has a terminal aldehyde group at the end of a long aliphatic chain,
    potentially including unsaturation (double bonds).

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
    # It looks for a carbonyl carbon (C=O) connected to a hydrogen atom (H)
    aldehyde_pattern = Chem.MolFromSmarts('[CX3H1](=O)')
    
    # Check the main features
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)

    if not aldehyde_matches:
        return False, "No terminal aldehyde group found"
    
    # Check for a sufficient aliphatic carbon chain attached to the aldehyde
    for match in aldehyde_matches:
        aldehyde_carbon_idx = match[0]  # The carbon atom in [CH]=O

        # Define criteria for valid aliphatic carbon atoms
        def is_valid_aliphatic(atom):
            if atom.GetAtomicNum() == 6 and not atom.IsInRing():  # Carbon atoms only, not in rings
                if atom.GetDegree() == 1 or (atom.GetDegree() == 2 and atom.GetIsAromatic() is False):
                    return True
            return False

        visited = set()
        chain_length = 0
        
        def traverse_chain(atom_idx, prev_idx):
            nonlocal chain_length
            if atom_idx not in visited:
                visited.add(atom_idx)
                atom = mol.GetAtomWithIdx(atom_idx)
                if is_valid_aliphatic(atom):
                    chain_length += 1
                    # Traverse neighbors except the one we came from
                    for neighbor in atom.GetNeighbors():
                        next_idx = neighbor.GetIdx()
                        if next_idx != prev_idx:
                            traverse_chain(next_idx, atom_idx)

        # Traverse from the aldehyde carbon towards the carbon chain
        carbon_neighbors = [nbr for nbr in mol.GetAtomWithIdx(aldehyde_carbon_idx).GetNeighbors() if nbr.GetAtomicNum() == 6]
        
        for neighbor in carbon_neighbors:
            traverse_chain(neighbor.GetIdx(), aldehyde_carbon_idx)

        if chain_length >= 6:  # Typical chain length criteria for fatty structures
            return True, "Valid fatty aldehyde: Terminal aldehyde group with a suitable aliphatic chain"
    
    return False, "Carbon chain too short for typical fatty aldehyde"