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

    # Define SMARTS pattern for N-sulfonylurea
    # [N;X3] - nitrogen with 3 connections
    # C(=O) - carbonyl group
    # [N;X3] - second nitrogen with 3 connections
    # S(=O)(=O) - sulfonyl group
    n_sulfonylurea_pattern = Chem.MolFromSmarts("[N;X3]C(=O)[N;X3]S(=O)(=O)")
    
    # Check for presence of N-sulfonylurea group
    if not mol.HasSubstructMatch(n_sulfonylurea_pattern):
        return False, "No N-sulfonylurea group found"
    
    # Additional checks to avoid false positives
    
    # Count number of N-sulfonylurea groups
    matches = mol.GetSubstructMatches(n_sulfonylurea_pattern)
    if len(matches) > 1:
        return True, f"Contains {len(matches)} N-sulfonylurea groups"
    
    # Verify essential atoms are present
    required_atoms = {'N', 'C', 'O', 'S'}
    mol_atoms = set(atom.GetSymbol() for atom in mol.GetAtoms())
    if not required_atoms.issubset(mol_atoms):
        return False, "Missing essential atoms for N-sulfonylurea"
    
    # Count sulfur atoms (should have at least one)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    if s_count == 0:
        return False, "No sulfur atoms found"
        
    # Count oxygen atoms (should have at least 3: one from urea, two from sulfonyl)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, "Insufficient oxygen atoms for N-sulfonylurea"
    
    # Additional pattern to ensure proper connectivity
    # This checks for the complete N-sulfonylurea group with surrounding atoms
    detailed_pattern = Chem.MolFromSmarts("[*]-[N;X3]C(=O)[N;X3]S(=O)(=O)-[*]")
    if not mol.HasSubstructMatch(detailed_pattern):
        return False, "N-sulfonylurea group lacks proper connectivity"

    return True, "Contains N-sulfonylurea group with proper connectivity"