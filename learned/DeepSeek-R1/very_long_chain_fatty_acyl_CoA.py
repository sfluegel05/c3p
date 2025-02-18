"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:XXXXX very long-chain fatty acyl-CoA (chain length > C22)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA (chain length > C22).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if very long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Match the CoA thioester pattern: C(=O)S-C-C-N...
    # Using SMARTS to find [C]=O connected to S, which is connected to C-C-N
    coa_thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2][CX4][CX4][NX2]")
    if not mol.HasSubstructMatch(coa_thioester_pattern):
        return False, "No CoA thioester group detected"
    
    # Find all thioester C=O-S groups
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found"
    
    max_chain_length = 0
    for match in thioester_matches:
        carbonyl_carbon = match[0]
        visited = set()
        stack = [(carbonyl_carbon, 0)]  # (atom index, current depth)
        current_max = 0
        
        while stack:
            atom_idx, depth = stack.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:  # Only consider carbon atoms in the chain
                continue
            
            current_max = max(current_max, depth + 1)  # depth starts at 0, so +1 for current atom
            # Traverse all neighboring carbons except those leading back to the thioester S
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                # Avoid going back towards the sulfur
                if neighbor.GetAtomicNum() == 16 and neighbor.GetIdx() == match[2]:
                    continue
                # Avoid oxygen in the carbonyl
                if neighbor.GetIdx() == match[1]:
                    continue
                stack.append((neighbor_idx, depth + 1))
        
        if current_max > max_chain_length:
            max_chain_length = current_max
    
    if max_chain_length > 22:
        return True, f"Fatty acyl chain length is {max_chain_length} (>C22)"
    else:
        return False, f"Fatty acyl chain length is {max_chain_length} (≤C22)"