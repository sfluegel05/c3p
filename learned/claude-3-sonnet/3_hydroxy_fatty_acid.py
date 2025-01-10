"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: 3-hydroxy fatty acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Specific pattern for 3-hydroxy fatty acids:
    # [OH]-[CH2/CH]-[CH2]-[C](=O)[OH]
    # This ensures the hydroxy group is specifically at the 3-position
    beta_hydroxy_pattern = Chem.MolFromSmarts("[OX2H1]-[CX4H1,CX4H2]-[CX4H2]-[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(beta_hydroxy_pattern):
        return False, "No hydroxy group at 3-position"

    # Exclude molecules with too many carboxylic acid groups
    if len(mol.GetSubstructMatches(carboxylic_pattern)) > 1:
        return False, "Multiple carboxylic acid groups present"

    # Count carbons in the main chain from the carboxylic acid
    def count_chain_carbons(mol, start_idx):
        visited = set()
        queue = [(start_idx, 0)]
        max_depth = 0
        
        while queue:
            current, depth = queue.pop(0)
            if current in visited:
                continue
            visited.add(current)
            max_depth = max(max_depth, depth)
            
            atom = mol.GetAtomWithIdx(current)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon
                    queue.append((neighbor.GetIdx(), depth + 1))
        
        return max_depth

    # Find the carboxylic carbon and count chain length
    matches = mol.GetSubstructMatches(carboxylic_pattern)
    if matches:
        carboxyl_carbon = matches[0][0]
        chain_length = count_chain_carbons(mol, carboxyl_carbon)
        if chain_length < 4:  # Minimum chain length for 3-hydroxy fatty acid
            return False, "Carbon chain too short"

    # Check for excessive branching
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            carbon_neighbors = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 6)
            if carbon_neighbors > 3:
                return False, "Excessive branching found"

    # Exclude cyclic structures except for allowed cyclopropyl groups
    ring_info = mol.GetRingInfo()
    for ring_size in ring_info.RingSizes():
        if ring_size != 3:  # Allow only cyclopropyl rings
            return False, "Contains non-cyclopropyl ring"

    # Count heteroatoms other than O
    heteroatoms = sum(1 for atom in mol.GetAtoms() 
                     if atom.GetAtomicNum() not in [1, 6, 8])
    if heteroatoms > 0:
        return False, "Contains non-allowed heteroatoms"

    # Check for aromatic systems
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings"

    # Count oxygens (should be 3-4 typically, but allow up to 6 for dihydroxy variants)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 3:
        return False, "Insufficient oxygen atoms"
    if oxygen_count > 6:
        return False, "Too many oxygen atoms"

    return True, "Contains 3-hydroxy group in fatty acid structure"