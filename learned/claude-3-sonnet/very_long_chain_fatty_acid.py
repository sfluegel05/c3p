"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: CHEBI:27208 very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    VLCFAs have chain length >C22. Those >C27 are ultra-long-chain fatty acids.

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

    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Check if molecule has large ring systems (likely not a fatty acid)
    ring_info = mol.GetRingInfo()
    if any(len(ring) > 7 for ring in ring_info.AtomRings()):
        return False, "Contains large ring systems - not a fatty acid"
    
    # Get the carboxyl carbon
    carboxyl_match = mol.GetSubstructMatch(carboxyl_pattern)
    carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_match[0])
    
    # Find the longest aliphatic chain starting from carboxyl group
    def find_aliphatic_chain(atom, visited=None, depth=0):
        if visited is None:
            visited = set()
            
        visited.add(atom.GetIdx())
        max_length = depth
        
        # Only consider carbon atoms that aren't in rings (unless small rings like cyclopropyl)
        for neighbor in atom.GetNeighbors():
            if (neighbor.GetIdx() not in visited and 
                neighbor.GetAtomicNum() == 6 and 
                (not neighbor.IsInRing() or 
                 (neighbor.IsInRing() and all(len(ring) <= 3 for ring in ring_info.AtomRings() if neighbor.GetIdx() in ring)))):
                length = find_aliphatic_chain(neighbor, visited.copy(), depth + 1)
                max_length = max(max_length, length)
                
        return max_length
    
    # Get the longest aliphatic chain length
    chain_length = find_aliphatic_chain(carboxyl_carbon)
    
    if chain_length <= 22:
        return False, f"Aliphatic chain length C{chain_length} is too short (must be >C22)"
    
    # Calculate molecular weight - fatty acids should be in a reasonable range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 1000:  # Most natural fatty acids are under 1000 Da
        return False, "Molecular weight too high for a fatty acid"
        
    # Check for fatty acid features (only count if connected to main chain)
    def is_connected_to_main_chain(atom):
        path_exists = False
        current = atom
        visited = set()
        while current and not path_exists:
            visited.add(current.GetIdx())
            if current.GetIdx() == carboxyl_match[0]:
                path_exists = True
                break
            for neighbor in current.GetNeighbors():
                if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum() == 6:
                    current = neighbor
                    break
            else:
                current = None
        return path_exists

    # Find features connected to main chain
    is_ultra = chain_length > 27
    
    double_bond_pattern = Chem.MolFromSmarts("[#6]=[#6]")
    has_double_bonds = any(is_connected_to_main_chain(mol.GetAtomWithIdx(match[0])) 
                          for match in mol.GetSubstructMatches(double_bond_pattern))
    
    triple_bond_pattern = Chem.MolFromSmarts("[#6]#[#6]")
    has_triple_bonds = any(is_connected_to_main_chain(mol.GetAtomWithIdx(match[0])) 
                          for match in mol.GetSubstructMatches(triple_bond_pattern))
    
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H1][#6]")
    has_hydroxy = any(is_connected_to_main_chain(mol.GetAtomWithIdx(match[1])) 
                      for match in mol.GetSubstructMatches(hydroxy_pattern))
    
    cyclopropyl_pattern = Chem.MolFromSmarts("[C]1[C][C]1")
    has_cyclopropyl = mol.HasSubstructMatch(cyclopropyl_pattern)
    
    branch_pattern = Chem.MolFromSmarts("[#6][CH]([#6])[#6]")
    has_branching = any(is_connected_to_main_chain(mol.GetAtomWithIdx(match[1])) 
                       for match in mol.GetSubstructMatches(branch_pattern))
    
    # Build description
    features = []
    if is_ultra:
        features.append("ultra-long-chain")
    if has_double_bonds:
        features.append("unsaturated")
    if has_triple_bonds:
        features.append("contains triple bonds")
    if has_hydroxy:
        features.append("hydroxylated")
    if has_cyclopropyl:
        features.append("contains cyclopropyl rings")
    if has_branching:
        features.append("branched")
    
    feature_str = " and ".join(features) if features else "saturated"
    chain_desc = f"C{chain_length}"
    
    return True, f"{chain_desc} {feature_str} fatty acid"