"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: CHEBI:27283 short-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def get_main_chain_length(mol, carboxyl_carbon_idx):
    """
    Get the length of the longest carbon chain starting from the carboxyl carbon.
    Uses a breadth-first search to find the longest continuous carbon chain.
    """
    visited = set([carboxyl_carbon_idx])
    current_level = {carboxyl_carbon_idx}
    chain_length = 1  # Start with 1 to count carboxyl carbon
    
    while current_level:
        next_level = set()
        for atom_idx in current_level:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                    next_level.add(neighbor.GetIdx())
                    visited.add(neighbor.GetIdx())
        if next_level:
            chain_length += 1
        current_level = next_level
    
    return chain_length

def has_invalid_substituents(mol):
    """
    Check if molecule has invalid substituents.
    Returns True if invalid substituents are found.
    Allows hydroxy (-OH) and oxo (=O) groups.
    """
    # Find carboxyl group atoms
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return True
    
    carboxyl_atoms = set(mol.GetSubstructMatch(carboxyl_pattern))
    
    # Allowed patterns (besides carboxyl)
    hydroxy_pattern = Chem.MolFromSmarts("[CX4][OX2H1]")  # -OH
    oxo_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")  # =O
    allowed_o_patterns = [hydroxy_pattern, oxo_pattern]
    
    # Get all oxygen atoms in allowed patterns
    allowed_o_atoms = set()
    for pattern in allowed_o_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                for atom_idx in match:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    if atom.GetAtomicNum() == 8:
                        allowed_o_atoms.add(atom_idx)
    
    # Check each atom
    for atom in mol.GetAtoms():
        # Skip hydrogens
        if atom.GetAtomicNum() == 1:
            continue
            
        # Only allow C, H, O atoms
        if atom.GetAtomicNum() not in [1, 6, 8]:
            return True
            
        # Check oxygen atoms
        if atom.GetAtomicNum() == 8:
            if atom.GetIdx() not in carboxyl_atoms and atom.GetIdx() not in allowed_o_atoms:
                return True
    
    return False

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    Short-chain fatty acids are aliphatic monocarboxylic acids with chain length < C6.
    
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
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) == 0:
        return False, "No carboxylic acid group found"
    if len(carboxyl_matches) > 1:
        return False, "Multiple carboxylic acid groups found"
        
    # Check for aromatic rings
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings"
        
    # Check for cyclic structures
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        return False, "Contains cyclic structures"
    
    # Check for invalid substituents
    if has_invalid_substituents(mol):
        return False, "Contains non-allowed substituents"
    
    # Get carbon chain length
    carboxyl_carbon_idx = carboxyl_matches[0][0]
    chain_length = get_main_chain_length(mol, carboxyl_carbon_idx)
    
    if chain_length < 2:
        return False, "Carbon chain too short"
    if chain_length > 5:
        return False, f"Carbon chain too long ({chain_length} carbons)"
        
    return True, f"Aliphatic monocarboxylic acid with {chain_length} carbons in longest chain"