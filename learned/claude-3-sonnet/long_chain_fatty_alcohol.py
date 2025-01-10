"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: long-chain fatty alcohol
Definition: A fatty alcohol with a chain length ranging from C13 to C22
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def count_chain_carbons(mol, start_idx, visited=None):
    """
    Recursively count connected carbons starting from a given atom
    """
    if visited is None:
        visited = set()
    
    if start_idx in visited:
        return 0
    
    visited.add(start_idx)
    atom = mol.GetAtomWithIdx(start_idx)
    
    if atom.GetAtomicNum() != 6:  # Not carbon
        return 0
        
    count = 1  # Count current carbon
    
    # Recursively explore neighbors
    for neighbor in atom.GetNeighbors():
        if neighbor.GetIdx() not in visited:
            count += count_chain_carbons(mol, neighbor.GetIdx(), visited)
            
    return count

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"
    
    # Get hydroxyl carbons (carbons attached to OH)
    hydroxyl_carbon_pattern = Chem.MolFromSmarts("[C][OH]")
    hydroxyl_carbons = mol.GetSubstructMatches(hydroxyl_carbon_pattern)
    
    if not hydroxyl_carbons:
        return False, "No carbon-hydroxyl bonds found"
    
    # Count total carbons
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if total_carbons < 13:
        return False, f"Too few total carbons (C{total_carbons}, need C13-C22)"
    if total_carbons > 22:
        # Allow slightly more carbons to account for branching
        if total_carbons > 25:
            return False, f"Too many carbons (C{total_carbons}, need C13-C22)"
    
    # Check for longest carbon chain from each hydroxyl carbon
    max_chain_length = 0
    for match in hydroxyl_carbons:
        chain_length = count_chain_carbons(mol, match[0])
        max_chain_length = max(max_chain_length, chain_length)
    
    if max_chain_length < 13:
        return False, f"Longest carbon chain too short (C{max_chain_length}, need C13-C22)"
    if max_chain_length > 22:
        return False, f"Longest carbon chain too long (C{max_chain_length}, need C13-C22)"
    
    # Check for excessive rings
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 2:  # Allow up to 2 rings
        return False, f"Too many rings ({ring_count}) for a fatty alcohol"
    
    # Check proportion of C and H atoms
    h_count = sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())
    total_atoms = mol.GetNumAtoms() + h_count
    ch_count = total_carbons + h_count
    
    if ch_count / total_atoms < 0.7:  # Allow more heteroatoms
        return False, "Not primarily a hydrocarbon structure"
    
    # Count other heteroatoms (excluding O from OH groups)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if n_count + s_count + p_count > 2:
        return False, "Too many heteroatoms for a fatty alcohol"
    
    return True, f"Contains hydroxyl group with appropriate carbon chain length (C{max_chain_length})"