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
    
    # Get the carbonyl and alpha carbon indices
    matches = mol.GetSubstructMatches(carboxylic_acid)
    carbonyl_idx = matches[0][0]
    alpha_idx = matches[0][1]
    
    # Traverse the hydrocarbon chain starting from alpha carbon, excluding the carbonyl
    visited = set()
    queue = [alpha_idx]
    branch_found = False
    
    while queue:
        current = queue.pop(0)
        if current in visited:
            continue
        visited.add(current)
        atom = mol.GetAtomWithIdx(current)
        
        # Count neighboring carbons excluding the carbonyl
        neighbor_carbons = []
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr.GetAtomicNum() == 6 and nbr_idx != carbonyl_idx:
                neighbor_carbons.append(nbr_idx)
        
        # If current atom has >=3 carbon neighbors, it's a branch point
        if len(neighbor_carbons) >= 3:
            branch_found = True
            break
        
        # Add neighboring carbons to the queue
        for nbr_idx in neighbor_carbons:
            if nbr_idx not in visited:
                queue.append(nbr_idx)
    
    if branch_found:
        return True, "Branched chain detected"
    else:
        return False, "No branched chain found"