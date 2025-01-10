"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
"""
Classifies: short-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for CoA structure components
    
    # Adenine base
    adenine_pattern = Chem.MolFromSmarts("c1nc(N)nc2c1ncn2")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine base (CoA component)"
    
    # Phosphate groups (at least 3)
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    phosphate_matches = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_matches < 3:
        return False, "Missing required phosphate groups (CoA component)"
        
    # Thioester linkage
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Missing thioester linkage"
        
    # Get the acyl part by splitting at the thioester
    # First, find all carbons connected to the thioester sulfur
    acyl_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:  # Sulfur
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon
                    acyl_carbons.append(neighbor.GetIdx())
                    
    if not acyl_carbons:
        return False, "Could not identify acyl chain"
        
    # Count carbons in the acyl chain
    # Create a visited set starting from the acyl carbon
    visited = set()
    stack = [acyl_carbons[0]]
    
    while stack:
        current = stack.pop()
        if current not in visited:
            visited.add(current)
            atom = mol.GetAtomWithIdx(current)
            if atom.GetAtomicNum() == 6:  # Carbon
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:  # Only follow carbon chain
                        stack.append(neighbor.GetIdx())
                        
    chain_length = len(visited)
    
    if chain_length > 6:
        return False, f"Acyl chain too long ({chain_length} carbons) for short-chain fatty acid"
    if chain_length < 2:
        return False, f"Acyl chain too short ({chain_length} carbons)"
        
    # Check for common modifications that are allowed
    allowed_modifications = [
        Chem.MolFromSmarts("[OX2H1]"), # hydroxy
        Chem.MolFromSmarts("[NX3H2]"), # amino
        Chem.MolFromSmarts("[CX3]=[CX3]"), # double bond
        Chem.MolFromSmarts("[CX4H3]") # methyl branch
    ]
    
    has_valid_modifications = True
    for mod in allowed_modifications:
        if mol.HasSubstructMatch(mod):
            has_valid_modifications = True
            break
            
    if not has_valid_modifications:
        return False, "No valid modifications found on acyl chain"
        
    return True, f"Valid short-chain ({chain_length} carbons) fatty acyl-CoA with appropriate modifications"