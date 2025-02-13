"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: CHEBI:27283 short-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def count_longest_carbon_chain(mol):
    """Helper function to count the longest carbon chain including the carboxyl carbon"""
    # Find carboxyl carbon first
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return 0
        
    # Get all carbons in molecule
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_indices:
        return 0
        
    # Get carboxyl carbon index
    carboxyl_match = mol.GetSubstructMatch(carboxyl_pattern)
    if not carboxyl_match:
        return 0
    carboxyl_carbon_idx = carboxyl_match[0]
    
    # Find longest path from carboxyl carbon
    max_length = 0
    for c_idx in carbon_indices:
        if c_idx != carboxyl_carbon_idx:
            path = Chem.GetShortestPath(mol, carboxyl_carbon_idx, c_idx)
            if path:
                # Count only carbons in path
                carbon_count = sum(1 for idx in path if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
                max_length = max(max_length, carbon_count)
    
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
        
    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Count number of carboxylic acid groups - should be exactly one
    carboxyl_matches = len(mol.GetSubstructMatches(carboxyl_pattern))
    if carboxyl_matches > 1:
        return False, "Multiple carboxylic acid groups found"
        
    # Check for aromatic rings - should have none
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings"
        
    # Check for non C,H,O atoms
    allowed_atoms = {6, 1, 8}  # C, H, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, f"Contains atoms other than C, H, O"
            
    # Get longest carbon chain length
    chain_length = count_longest_carbon_chain(mol)
    if chain_length < 2:
        return False, "Carbon chain too short"
    if chain_length > 5:
        return False, f"Carbon chain too long ({chain_length} carbons)"
        
    # Check for cyclic structures
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        return False, "Contains cyclic structures"
    
    return True, f"Aliphatic monocarboxylic acid with {chain_length} carbons in longest chain"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:27283',
        'name': 'short-chain fatty acid',
        'definition': 'An aliphatic monocarboxylic acid with a chain length of less than C6. '
                     'If any non-hydrocarbon substituent is present, the compound is not normally '
                     'regarded as a short-chain fatty acid.',
    }
}