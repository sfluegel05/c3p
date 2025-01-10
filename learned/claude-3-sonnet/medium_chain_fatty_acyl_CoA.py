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
    Medium-chain fatty acids typically have 4-12 carbons in a primarily linear chain.
    
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
    adenine_pattern = Chem.MolFromSmarts("[n,N]1[c,C][n,N][c,C]2[c,C]([N,n])[n,N][c,C][n,N][c,C]12")
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)O")
    pantetheine_pattern = Chem.MolFromSmarts("CCNC(=O)CCNC(=O)C(O)C(C)(C)COP")
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    
    # Basic structure checks
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine moiety"
    
    phosphate_matches = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_matches < 2:
        return False, f"Insufficient phosphate groups (found {phosphate_matches}, need at least 2)"
    
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing pantetheine moiety"
    
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"
    
    # Get the carbon atom index of the thioester carbonyl
    carbonyl_idx = thioester_matches[0][0]
    
    # Patterns that should not be present in fatty acids
    aromatic_pattern = Chem.MolFromSmarts("a")
    if mol.HasSubstructMatch(aromatic_pattern):
        return False, "Contains aromatic rings (not a fatty acid)"
    
    # Count carbons in the chain using modified BFS that follows linear chains
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
                # Only count carbons that are part of the main chain or single methyl branches
                if len(list(current_atom.GetNeighbors())) <= 3:  # Allow branching
                    chain_carbons.add(current_idx)
                
            # Add neighbors to queue if they're carbon and not part of a ring
            for neighbor in current_atom.GetNeighbors():
                if (neighbor.GetAtomicNum() == 6  # is carbon
                    and neighbor.GetIdx() not in visited
                    and not neighbor.IsInRing()  # exclude ring carbons
                    and len(list(neighbor.GetNeighbors())) <= 3):  # max 3 neighbors for linear chain
                    queue.append(neighbor.GetIdx())
        
        return len(chain_carbons)
    
    chain_length = count_chain_carbons(mol, carbonyl_idx)
    
    # Check chain length (4-12 carbons for medium chain)
    if chain_length < 4:
        return False, f"Fatty acid chain too short ({chain_length} carbons)"
    if chain_length > 12:
        return False, f"Fatty acid chain too long ({chain_length} carbons)"
    
    # Additional checks for reasonable molecular properties
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 700 or mol_weight > 1200:
        return False, f"Molecular weight {mol_weight:.1f} outside typical range"
    
    # Check for excessive ring systems
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 2:  # Allow for adenine ring system
        return False, "Contains too many ring systems"
        
    return True, f"Medium-chain fatty acyl-CoA with {chain_length} carbons in fatty acid chain"