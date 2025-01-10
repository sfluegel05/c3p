"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: methyl-branched fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    A methyl-branched fatty acid is a fatty acid containing only methyl branches.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
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
    
    # Count carboxylic acid groups - should have just one
    carboxylic_matches = len(mol.GetSubstructMatches(carboxylic_pattern))
    if carboxylic_matches != 1:
        return False, f"Invalid number of carboxylic acid groups ({carboxylic_matches})"
    
    # Check for cyclic structures (excluding small rings that might be part of side groups)
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(r) for r in ring_info.AtomRings()]
    if any(size >= 6 for size in ring_sizes):
        return False, "Contains large ring structure"
    
    # Count carbons and check if it's a reasonable size for a fatty acid
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Carbon chain too short"
    if c_count > 30:
        return False, "Carbon chain too long for typical fatty acid"
        
    # Check for aromatic character
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings"
    
    # Count heteroatoms (excluding O from COOH)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if s_count > 0 or p_count > 0:
        return False, "Contains unexpected heteroatoms (S, P)"
    if n_count > 1:  # More restrictive on N
        return False, "Too many nitrogen atoms"
    if o_count > 3:  # Allow COOH (2 O) plus maybe one extra O
        return False, "Too many oxygen atoms"
    
    # Look for methyl branches using multiple patterns
    methyl_patterns = [
        "[CH3][CX4]", # Standard methyl branch
        "[CH3][CH1]", # Methyl on branching point
        "[CH3][C]([CH3])[CH3]" # Geminal dimethyl
    ]
    
    total_methyl_branches = 0
    for pattern in methyl_patterns:
        matches = len(mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)))
        total_methyl_branches += matches
    
    if total_methyl_branches == 0:
        return False, "No methyl branches found"
    
    # Check for non-methyl branches (e.g. ethyl, propyl)
    # This pattern matches carbon chains longer than methyl attached to the main chain
    non_methyl_branch = Chem.MolFromSmarts("[CH2][CH2][CH2,CH3]")
    if mol.HasSubstructMatch(non_methyl_branch):
        # Verify if it's part of the main chain rather than a branch
        matches = len(mol.GetSubstructMatches(non_methyl_branch))
        if matches > 1:  # Allow one match which might be the main chain
            return False, "Contains non-methyl branches"
    
    # Success case
    return True, f"Contains carboxylic acid group with methyl branching"