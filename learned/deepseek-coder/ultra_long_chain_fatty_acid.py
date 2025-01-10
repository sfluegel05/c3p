"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
"""
Classifies: CHEBI:XXXXX ultra-long-chain fatty acid
"""
from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.
    An ultra-long-chain fatty acid is defined as having a carbon chain length greater than C27.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ultra-long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"

    # Find the longest carbon chain
    longest_chain_length = 0
    best_chain = []
    
    # Iterate through all carboxylic acid groups
    for match in carboxylic_acid_matches:
        # Start from the carbon in the carboxylic acid group
        start_atom = mol.GetAtomWithIdx(match[0])
        visited = set()
        stack = [(start_atom, 1, [start_atom.GetIdx()])]  # (current_atom, current_length, current_chain)
        
        while stack:
            current_atom, current_length, current_chain = stack.pop()
            visited.add(current_atom.GetIdx())
            
            # Update longest chain if necessary
            if current_length > longest_chain_length:
                longest_chain_length = current_length
                best_chain = current_chain
            
            # Explore neighbors
            for neighbor in current_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                    # Continue if the neighbor is not already in the chain
                    if neighbor.GetIdx() not in current_chain:
                        stack.append((neighbor, current_length + 1, current_chain + [neighbor.GetIdx()]))

    # Check if the longest chain is greater than C27
    if longest_chain_length > 27:
        return True, f"Longest carbon chain length is {longest_chain_length} (greater than C27)"
    else:
        return False, f"Longest carbon chain length is {longest_chain_length} (not greater than C27)"