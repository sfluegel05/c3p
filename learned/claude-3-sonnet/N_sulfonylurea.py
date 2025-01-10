"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
"""
Classifies: N-sulfonylurea compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    N-sulfonylureas contain a urea group where one nitrogen is connected to a sulfonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-sulfonylurea, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for N-sulfonylurea
    # Pattern 1: Sulfonyl connected to first nitrogen
    pattern1 = Chem.MolFromSmarts("[N]S(=O)(=O)[*]")
    # Pattern 2: Urea group
    pattern2 = Chem.MolFromSmarts("[N]C(=O)[N]")
    
    if not (mol.HasSubstructMatch(pattern1) and mol.HasSubstructMatch(pattern2)):
        return False, "Missing required N-sulfonylurea substructure"

    # More specific pattern to ensure proper connectivity
    # Allows for both possible arrangements and various nitrogen states
    nsulf_pattern1 = Chem.MolFromSmarts("[#7]S(=O)(=O)[*]C(=O)[#7]")
    nsulf_pattern2 = Chem.MolFromSmarts("[#7]C(=O)[#7]S(=O)(=O)[*]")
    
    has_arrangement1 = mol.HasSubstructMatch(nsulf_pattern1)
    has_arrangement2 = mol.HasSubstructMatch(nsulf_pattern2)
    
    if not (has_arrangement1 or has_arrangement2):
        return False, "N-sulfonylurea group lacks proper connectivity"

    # Additional validation to avoid false positives
    
    # Check that carbons and nitrogens are properly connected
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            # Check nitrogen's environment
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if 'S' in neighbors and 'C' in neighbors:
                # Verify carbon is part of carbonyl group
                for n in atom.GetNeighbors():
                    if n.GetSymbol() == 'C':
                        c_neighbors = [nb.GetSymbol() for nb in n.GetNeighbors()]
                        if 'O' in c_neighbors and c_neighbors.count('N') >= 1:
                            return True, "Contains N-sulfonylurea group with validated connectivity"

    # Count key atoms to ensure proper composition
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if s_count < 1 or o_count < 3 or n_count < 2:
        return False, "Insufficient atoms for N-sulfonylurea structure"

    # Pattern to exclude similar but invalid structures
    invalid_pattern = Chem.MolFromSmarts("[N]S(=O)(=O)[N]")
    if mol.HasSubstructMatch(invalid_pattern):
        # Check if this is a direct S-N-N connection (invalid) vs proper N-sulfonylurea
        for match in mol.GetSubstructMatches(invalid_pattern):
            n1, s, n2 = match
            if len([n for n in mol.GetAtomWithIdx(n1).GetNeighbors() 
                   if n.GetSymbol() == 'C' and any(nb.GetSymbol() == 'O' 
                   for nb in n.GetNeighbors())]) == 0:
                return False, "Contains invalid sulfonamide arrangement"

    return True, "Contains N-sulfonylurea group with proper connectivity"