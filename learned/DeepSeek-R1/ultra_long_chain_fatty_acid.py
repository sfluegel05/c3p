"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
"""
Classifies: CHEBI: ??? ultra-long-chain fatty acid (chain length > C27)
"""
from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid (chain length > C27).
    Must contain a carboxylic acid group and have a main carbon chain longer than 27 carbons.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Match carboxylic acid group (-COOH)
    carboxylic_acid = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group"
    
    # Get all carboxylic acid matches
    matches = mol.GetSubstructMatches(carboxylic_acid)
    
    # Helper to calculate longest carbon chain from start atom
    def get_chain_length(start_atom):
        visited = set()
        stack = [(start_atom, 1)]  # (atom, current_length)
        max_len = 0
        
        while stack:
            atom, length = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            max_len = max(max_len, length)
            
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                    stack.append((neighbor, length + 1))
        
        return max_len
    
    # Check each carboxylic acid group's chain length
    for match in matches:
        carbonyl_c = mol.GetAtomWithIdx(match[0])
        # Find alpha carbon (non-oxygen neighbor)
        alpha_carbons = [n for n in carbonyl_c.GetNeighbors() if n.GetAtomicNum() == 6]
        if not alpha_carbons:
            continue
        
        chain_length = get_chain_length(alpha_carbons[0])
        # Total chain length = alpha chain + carbonyl carbon
        if chain_length >= 27:
            return True, f"Main chain length {chain_length + 1} > C27"
    
    return False, "Main chain length ≤ C27"