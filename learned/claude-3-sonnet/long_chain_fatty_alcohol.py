"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: long-chain fatty alcohol
Definition: A fatty alcohol with a chain length ranging from C13 to C22
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.GraphDescriptors import CalcNumRings

def find_longest_chain(mol, start_idx, exclude_atoms=None):
    """
    Find the longest carbon chain from a starting atom using DFS
    """
    if exclude_atoms is None:
        exclude_atoms = set()
    
    def dfs(current_idx, path):
        nonlocal longest_path
        
        if len(path) > len(longest_path):
            longest_path = path[:]
            
        atom = mol.GetAtomWithIdx(current_idx)
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if (n_idx not in path and 
                n_idx not in exclude_atoms and 
                neighbor.GetAtomicNum() == 6):  # Carbon atoms only
                dfs(n_idx, path + [n_idx])
    
    longest_path = []
    dfs(start_idx, [start_idx])
    return len(longest_path)

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

    # Exclude carboxylic acids
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if mol.HasSubstructMatch(carboxyl_pattern):
        return False, "Contains carboxylic acid group"
    
    # Find primary and secondary alcohols
    primary_alcohol = Chem.MolFromSmarts("[CH2][OH]")
    secondary_alcohol = Chem.MolFromSmarts("[CH1][OH]")
    
    primary_matches = mol.GetSubstructMatches(primary_alcohol)
    secondary_matches = mol.GetSubstructMatches(secondary_alcohol)
    
    if not (primary_matches or secondary_matches):
        return False, "No primary or secondary alcohol group found"
    
    # Count rings - fatty alcohols should be mainly aliphatic
    if CalcNumRings(mol) > 1:
        return False, "Too many rings for a fatty alcohol"
    
    # Get all OH-attached carbons
    all_alcohol_matches = list(primary_matches) + list(secondary_matches)
    
    # Find longest carbon chain from each OH-attached carbon
    max_chain_length = 0
    for match in all_alcohol_matches:
        c_idx = match[0]  # Carbon attached to OH
        chain_length = find_longest_chain(mol, c_idx)
        max_chain_length = max(max_chain_length, chain_length)
    
    if max_chain_length < 13:
        return False, f"Carbon chain too short (C{max_chain_length}, need C13-C22)"
    if max_chain_length > 22:
        return False, f"Carbon chain too long (C{max_chain_length}, need C13-C22)"
    
    # Check for excessive heteroatoms
    n_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[#7]")))
    s_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[#16]")))
    p_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[#15]")))
    
    if n_count + s_count + p_count > 1:
        return False, "Too many heteroatoms for a fatty alcohol"
    
    # Count carbons and hydrogens to ensure mainly hydrocarbon nature
    c_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[#6]")))
    if c_count < 13:
        return False, f"Insufficient carbons (C{c_count}, need C13-C22)"
        
    return True, f"Contains alcohol group with appropriate carbon chain length (C{max_chain_length})"