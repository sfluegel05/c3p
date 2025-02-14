"""
Classifies: CHEBI:24026 fatty alcohol
"""
"""
Classifies: fatty alcohol
"""
from rdkit import Chem

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is an aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms.
    Fatty alcohols may be saturated or unsaturated and may be branched or unbranched.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify hydroxyl groups (OH) attached to carbon atoms
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4,CX3H1,CX2H2,CX1H3]-[OX2H]")
    matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    if not matches:
        return False, "No hydroxyl group attached to carbon found"
    
    # Check each hydroxyl group
    for match in matches:
        carbon_idx = match[0]
        oxygen_idx = match[1]
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        
        # Perform DFS to count connected aliphatic carbons
        visited = set()
        stack = [carbon_idx]
        carbon_count = 0

        while stack:
            atom_idx = stack.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            
            # Consider only aliphatic carbons
            if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic():
                carbon_count += 1
                # Add neighboring atoms to stack
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited and neighbor.GetAtomicNum() == 6 and not neighbor.GetIsAromatic():
                        stack.append(neighbor_idx)

        # Check if the carbon chain length is at least 3
        if carbon_count >= 3:
            return True, f"Contains aliphatic alcohol group with a chain of {carbon_count} carbons"

    return False, "No aliphatic chain with at least 3 carbons attached to hydroxyl group"

# Example usage:
# smiles = "CCCCCCCCCCCCCC(O)CCC"  # tetradecan-5-ol
# result, reason = is_fatty_alcohol(smiles)
# print(result, reason)