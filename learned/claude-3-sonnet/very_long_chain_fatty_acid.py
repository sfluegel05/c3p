"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdmolops import GetShortestPath

def count_chain_carbons(mol):
    """
    Counts the effective chain length by finding the longest carbon chain
    that includes the carboxylic acid carbon.
    """
    # Find the carboxylic acid carbon
    carboxyl_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return 0
    
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    carboxyl_carbon = carboxyl_matches[0][0]
    
    # Get all carbon atoms
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() 
                     if atom.GetAtomicNum() == 6]
    
    # Find the longest path from carboxyl carbon to any other carbon
    max_length = 0
    for carbon_idx in carbon_indices:
        path_length = len(GetShortestPath(mol, carboxyl_carbon, carbon_idx))
        max_length = max(max_length, path_length)
    
    return max_length

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid (>C22).
    
    Args:
        smiles: SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) indicating if molecule is VLCFA and reason
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carbons in longest chain
    chain_length = count_chain_carbons(mol)
    if chain_length < 22:
        return False, f"Chain length ({chain_length}) less than C22"
    
    # Basic element counts
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    # Check for reasonable element ratios
    if o_count > 8:  # Allow for some hydroxy groups
        return False, "Too many oxygen atoms for a fatty acid"
    if n_count > 2:  # Allow for some N-containing modifications
        return False, "Too many nitrogen atoms for a fatty acid"
        
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 700 and n_count == 0:  # Higher limit for N-containing compounds
        return False, "Molecular weight too high for typical VLCFA"
        
    # Count rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 3:
        return False, "Too many rings for a fatty acid"
    
    # Check for sugar-like patterns (too many OH groups in close proximity)
    polyol_pattern = Chem.MolFromSmarts('[OH][CH][CH][OH]')
    if mol.HasSubstructMatch(polyol_pattern):
        return False, "Contains sugar-like polyol pattern"
        
    # Classify as ultra-long-chain if >C27
    if chain_length > 27:
        return True, f"Ultra-long-chain fatty acid (C{chain_length})"
    else:
        return True, f"Very long-chain fatty acid (C{chain_length})"