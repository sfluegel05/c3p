"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: CHEBI:90255 epoxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid contains a carboxylic acid group and an epoxide ring (3-membered cyclic ether).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for carboxylic acid group (COOH or COO-)
    # Using more flexible pattern to account for different protonation states
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1,O-]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group"
    
    # Check for epoxide ring (3-membered cyclic ether: C-O-C)
    # Corrected SMARTS pattern for 3-membered epoxide ring
    epoxide_pattern = Chem.MolFromSmarts("[O]1CC1")
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide ring found"
    
    # Validate fatty acid characteristics
    # Check chain length: at least 12 carbons in the longest chain
    # More accurate than total carbon count
    chain_length = 0
    longest_chain = 0
    
    # Get all chains starting from the carboxylic acid
    for match in mol.GetSubstructMatches(carboxyl_pattern):
        acid_carbon = match[0]  # Carbon in CX3(=O)
        # Traverse along the carbon chain
        current_chain = 1  # Start counting from acid carbon
        visited = set()
        stack = [(acid_carbon, current_chain)]
        
        while stack:
            atom_idx, count = stack.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:  # Carbon atom
                for neighbor in atom.GetNeighbors():
                    # Avoid backtracking through oxygen atoms (ester linkages etc.)
                    if neighbor.GetAtomicNum() == 8:
                        continue
                    if neighbor.GetIdx() not in visited:
                        stack.append((neighbor.GetIdx(), count + 1))
                        if count + 1 > longest_chain:
                            longest_chain = count + 1
    
    if longest_chain < 12:
        return False, f"Main chain too short ({longest_chain} carbons), need ≥12"
    
    return True, "Contains carboxylic acid and 3-membered epoxide ring with sufficient chain length"