"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: medium-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.
    Medium-chain fatty acids typically have 6-12 carbons.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for CoA moiety patterns
    # More flexible adenine pattern matching both forms seen in examples
    adenine_pattern = Chem.MolFromSmarts("[n,N]1[c,C][n,N][c,C]2[c,C]([N,n])[n,N][c,C][n,N][c,C]12")
    
    # Phosphate groups pattern
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)O")
    
    # Pantetheine pattern with thiol
    pantetheine_pattern = Chem.MolFromSmarts("CCNC(=O)CCNC(=O)C(O)C(C)(C)COP")
    
    # Thioester linkage pattern
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    
    # Check for key structural elements
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine moiety"
    
    phosphate_matches = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_matches < 2:
        return False, f"Insufficient phosphate groups (found {phosphate_matches}, need at least 2)"
    
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing pantetheine moiety"
    
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"
    
    # Find the fatty acid chain
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Could not analyze fatty acid chain"
    
    # Get the carbon atom index of the thioester carbonyl
    carbonyl_idx = thioester_matches[0][0]
    
    # Count carbons in the chain using BFS
    def count_chain_carbons(mol, start_idx):
        visited = set()
        chain_carbons = set()
        queue = [start_idx]
        
        while queue:
            current_idx = queue.pop(0)
            if current_idx in visited:
                continue
                
            visited.add(current_idx)
            current_atom = mol.GetAtomWithIdx(current_idx)
            
            # Count carbons (excluding the carbonyl carbon)
            if current_atom.GetAtomicNum() == 6 and current_idx != start_idx:
                chain_carbons.add(current_idx)
                
            # Add neighbors to queue if they're not sulfur and not oxygen
            for neighbor in current_atom.GetNeighbors():
                if (neighbor.GetAtomicNum() not in [8, 16]  # not oxygen or sulfur
                    and neighbor.GetIdx() not in visited):
                    queue.append(neighbor.GetIdx())
        
        return len(chain_carbons)
    
    chain_length = count_chain_carbons(mol, carbonyl_idx)
    
    # Check chain length (6-12 carbons for medium chain)
    if chain_length < 6:
        return False, f"Fatty acid chain too short ({chain_length} carbons)"
    if chain_length > 12:
        return False, f"Fatty acid chain too long ({chain_length} carbons)"
    
    # Additional checks for reasonable molecular properties
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 700 or mol_weight > 1200:
        return False, f"Molecular weight {mol_weight:.1f} outside typical range"
    
    return True, f"Medium-chain fatty acyl-CoA with {chain_length} carbons in fatty acid chain"