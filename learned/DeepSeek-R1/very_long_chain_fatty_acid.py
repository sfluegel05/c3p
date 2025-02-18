"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid (chain length > C22) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a VLCFA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find carboxylic acid group (-COOH)
    carb_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carb_acid_pattern):
        return False, "No carboxylic acid group found"

    # Get matches for carboxylic acid
    matches = mol.GetSubstructMatches(carb_acid_pattern)
    max_chain_length = 0

    for match in matches:
        carbonyl_carbon = match[0]
        # Find adjacent carbon (alpha carbon) to start chain measurement
        alpha_carbons = [n for n in mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors() 
                        if n.GetAtomicNum() == 6 and n.GetIdx() != match[1]]
        
        if not alpha_carbons:
            continue
            
        # Perform BFS to find longest carbon chain from alpha carbon
        visited = set()
        queue = [(alpha_carbons[0], 1)]  # (atom object, chain length)
        current_max = 0
        
        while queue:
            atom, length = queue.pop(0)
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            
            if length > current_max:
                current_max = length
                
            # Add neighboring carbons (allow single/double/triple bonds)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                    queue.append((neighbor, length + 1))

        # Total chain length includes the carbonyl carbon (+1)
        total_length = current_max + 1
        if total_length > max_chain_length:
            max_chain_length = total_length

    if max_chain_length > 22:
        return True, f"Main carbon chain length ({max_chain_length}) > 22"
    else:
        return False, f"Main carbon chain length ({max_chain_length}) ≤ 22"