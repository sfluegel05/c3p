"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
"""
Classifies: CHEBI_59820 branched-chain fatty acid
"""
from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    A branched-chain fatty acid has a hydrocarbon chain with one or more alkyl substituents.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group"
    
    # Get carbonyl and alpha carbon indices
    matches = mol.GetSubstructMatches(carboxylic_acid)
    carbonyl_idx = matches[0][0]
    alpha_idx = matches[0][1]
    
    # First check alpha carbon for immediate branching
    alpha_atom = mol.GetAtomWithIdx(alpha_idx)
    alpha_neighbors = [nbr.GetIdx() for nbr in alpha_atom.GetNeighbors() 
                      if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != carbonyl_idx]
    
    # If alpha carbon has more than 1 carbon neighbor (excluding carbonyl), it's branched
    if len(alpha_neighbors) >= 2:
        return True, "Branch at alpha carbon"
    
    # Then check rest of the chain for branching points
    visited = set([alpha_idx, carbonyl_idx])
    queue = list(alpha_neighbors)  # Start with alpha's neighbors
    
    while queue:
        current_idx = queue.pop(0)
        if current_idx in visited:
            continue
        visited.add(current_idx)
        
        current_atom = mol.GetAtomWithIdx(current_idx)
        neighbors = [nbr.GetIdx() for nbr in current_atom.GetNeighbors() 
                    if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != carbonyl_idx]
        
        # Branching occurs if a carbon has 3+ connections (including chain continuation)
        if len(neighbors) >= 3:
            return True, "Branch in hydrocarbon chain"
        
        # Add unvisited carbon neighbors to queue
        for nbr_idx in neighbors:
            if nbr_idx not in visited:
                queue.append(nbr_idx)
    
    return False, "No branched chain found"