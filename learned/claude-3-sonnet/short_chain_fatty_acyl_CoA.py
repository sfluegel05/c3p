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
    
    # Adenine base - using more flexible SMARTS patterns
    adenine_patterns = [
        Chem.MolFromSmarts("n1cnc2c(N)ncnc12"),  # Most common form
        Chem.MolFromSmarts("N1C=NC2=C1N=CN=C2N"),  # Alternative form
        Chem.MolFromSmarts("[#7]1:[#6]:[#7]:[#6]2:[#6]1:[#7]:[#6]:[#7]:[#6]:2[#7]")  # Generic form
    ]
    has_adenine = False
    for pattern in adenine_patterns:
        if mol.HasSubstructMatch(pattern):
            has_adenine = True
            break
    if not has_adenine:
        return False, "Missing adenine base (CoA component)"
    
    # Pantetheine arm with thioester
    pantetheine_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)CCSC(=O)")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing pantetheine arm with thioester"
    
    # Phosphate groups (at least 3)
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    phosphate_matches = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_matches < 3:
        return False, f"Insufficient phosphate groups (found {phosphate_matches}, need ≥3)"
        
    # Get the acyl part by splitting at the thioester
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Missing thioester linkage"
    
    # Find the acyl carbon (the one double bonded to oxygen in thioester)
    acyl_carbon = thioester_matches[0][0]
    
    # Count carbons in the acyl chain using BFS
    visited = set()
    stack = [acyl_carbon]
    while stack:
        current = stack.pop()
        if current not in visited:
            visited.add(current)
            atom = mol.GetAtomWithIdx(current)
            if atom.GetAtomicNum() == 6:  # Carbon
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                        stack.append(neighbor.GetIdx())
    
    chain_length = len([x for x in visited if mol.GetAtomWithIdx(x).GetAtomicNum() == 6])
    
    if chain_length > 6:
        return False, f"Acyl chain too long ({chain_length} carbons) for short-chain fatty acid"
    if chain_length < 2:
        return False, f"Acyl chain too short ({chain_length} carbons)"
    
    # Check for common modifications
    modifications = []
    mod_patterns = {
        "hydroxy": ("[OX2H1]", "hydroxyl"),
        "amino": ("[NX3H2]", "amino"),
        "double_bond": ("[CX3]=[CX3]", "unsaturated"),
        "methyl_branch": ("[CX4H3]", "methyl branch")
    }
    
    for substructure, (pattern, name) in mod_patterns.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            modifications.append(name)
    
    mod_text = " with " + ", ".join(modifications) if modifications else ""
    return True, f"Valid short-chain ({chain_length} carbons) fatty acyl-CoA{mod_text}"