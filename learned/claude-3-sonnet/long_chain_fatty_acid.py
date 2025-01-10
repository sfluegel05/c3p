"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: Long-chain fatty acids (C13-C22)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.GraphDescriptors import LongestPath

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid (C13-C22) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts('C(=O)[OH]')
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Get all carboxylic acid carbons
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) > 1:
        return False, "Multiple carboxylic acid groups found"
    
    # Find the carboxylic acid carbon
    carboxyl_carbon = carboxyl_matches[0][0]
    
    def count_carbon_chain(mol, start_atom, visited=None):
        """Helper function to count longest carbon chain from starting atom"""
        if visited is None:
            visited = set()
            
        visited.add(start_atom)
        max_length = 0
        
        for neighbor in start_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:  # Carbon atoms only
                length = count_carbon_chain(mol, neighbor, visited.copy())
                max_length = max(max_length, length)
                
        return max_length + 1
    
    # Count carbons in longest chain containing carboxylic group
    start_atom = mol.GetAtomWithIdx(carboxyl_carbon)
    chain_length = count_carbon_chain(mol, start_atom)
    
    # Check chain length (13-22 carbons)
    if chain_length < 13:
        return False, f"Carbon chain too short ({chain_length} carbons, need 13-22)"
    if chain_length > 22:
        return False, f"Carbon chain too long ({chain_length} carbons, need 13-22)"
        
    # Additional checks for reasonable molecular composition
    # Count total carbons and oxygens
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Most fatty acids shouldn't have more than 4 oxygens unless they're heavily modified
    if num_oxygens > 6:
        return False, "Too many oxygen atoms for typical fatty acid"
    
    # The main chain should constitute the majority of carbons
    if num_carbons > chain_length * 2:
        return False, "Too many carbons outside main chain"
        
    return True, f"Long-chain fatty acid with {chain_length} carbons in main chain"