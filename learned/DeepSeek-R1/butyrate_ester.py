"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: CHEBI:77518 butyrate ester
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester has a butyryl group (CH2CH2CH2CO-) attached via an ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all ester groups (O connected to carbonyl and another carbon)
    ester_pattern = Chem.MolFromSmarts("[OX2][C](=O)[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    for match in ester_matches:
        # Get the carbonyl carbon index (second atom in the match)
        carbonyl_idx = match[1]
        
        # Traverse the acyl chain starting from carbonyl carbon, excluding the ester oxygen
        visited = set()
        stack = [carbonyl_idx]
        carbon_count = 0
        
        while stack:
            current_idx = stack.pop()
            if current_idx in visited:
                continue
            visited.add(current_idx)
            current_atom = mol.GetAtomWithIdx(current_idx)
            
            # Only count carbons
            if current_atom.GetAtomicNum() != 6:
                continue
            
            carbon_count += 1
            
            # Add neighboring carbons except the ester oxygen
            for neighbor in current_atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                # Skip the ester oxygen (match[0])
                if neighbor_idx == match[0]:
                    continue
                if neighbor_idx not in visited:
                    stack.append(neighbor_idx)
        
        # Subtract 1 to exclude the carbonyl carbon itself (we want carbons in the acyl chain)
        acyl_chain_length = carbon_count - 1
        
        # Butyrate has 3 carbons in the acyl chain (CH2CH2CH2-)
        if acyl_chain_length == 3:
            return True, "Contains butyryl group (CCCC=O) connected via ester bond"
    
    return False, "No butyrate ester group detected"