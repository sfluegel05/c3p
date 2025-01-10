"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: CHEBI:27283 short-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def has_non_hydrocarbon_substituents(mol):
    """Check if molecule has any non-hydrocarbon substituents besides the carboxyl group"""
    # Find carboxyl group atoms
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return True
    
    carboxyl_atoms = set(mol.GetSubstructMatch(carboxyl_pattern))
    
    # Check each atom
    for atom in mol.GetAtoms():
        # Skip hydrogens and atoms in carboxyl group
        if atom.GetAtomicNum() == 1 or atom.GetIdx() in carboxyl_atoms:
            continue
        
        # Only allow carbon atoms
        if atom.GetAtomicNum() != 6:
            return True
            
        # Check for any oxygen attachments (except in carboxyl group)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() not in carboxyl_atoms:
                return True
    
    return False

def get_carbon_chain_length(mol):
    """Get the length of the main carbon chain including the carboxyl carbon"""
    # Find carboxyl carbon
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return 0
    
    carboxyl_match = mol.GetSubstructMatch(carboxyl_pattern)
    carboxyl_carbon = carboxyl_match[0]
    
    # Get all carbons
    carbon_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    
    # Find longest chain from carboxyl carbon
    max_length = 1  # Start with 1 to count carboxyl carbon
    visited = set([carboxyl_carbon])
    
    def dfs(atom_idx, current_length):
        nonlocal max_length
        max_length = max(max_length, current_length)
        
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                visited.add(neighbor.GetIdx())
                dfs(neighbor.GetIdx(), current_length + 1)
                visited.remove(neighbor.GetIdx())
    
    dfs(carboxyl_carbon, 1)
    return max_length

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for exactly one carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    carboxyl_matches = len(mol.GetSubstructMatches(carboxyl_pattern))
    if carboxyl_matches == 0:
        return False, "No carboxylic acid group found"
    if carboxyl_matches > 1:
        return False, "Multiple carboxylic acid groups found"
        
    # Check for aromatic rings
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings"
        
    # Check for non C,H,O atoms
    allowed_atoms = {6, 1, 8}  # C, H, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, f"Contains atoms other than C, H, O"
            
    # Check for cyclic structures
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        return False, "Contains cyclic structures"
    
    # Check for non-hydrocarbon substituents
    if has_non_hydrocarbon_substituents(mol):
        return False, "Contains non-hydrocarbon substituents"
    
    # Get carbon chain length
    chain_length = get_carbon_chain_length(mol)
    if chain_length < 2:
        return False, "Carbon chain too short"
    if chain_length >= 6:
        return False, f"Carbon chain too long ({chain_length} carbons)"
        
    return True, f"Aliphatic monocarboxylic acid with {chain_length} carbons in longest chain"