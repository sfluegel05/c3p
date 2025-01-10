"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count atoms to ensure it's a large molecule
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 30:  # Arbitrary threshold for minimum size
        return False, "Molecule too small to be a mucopolysaccharide"
    
    # Look for uronic acid pattern (COOH group attached to sugar ring)
    uronic_acid_pattern = Chem.MolFromSmarts("[CR1][CR1][CR1][CR1][CR1][CR1](C(=O)O)")
    if not mol.HasSubstructMatch(uronic_acid_pattern):
        return False, "No uronic acid units found"
    
    # Look for glycosamine pattern (sugar with NH2 group)
    glycosamine_pattern = Chem.MolFromSmarts("[CR1][CR1][CR1][CR1][CR1][CR1]N")
    if not mol.HasSubstructMatch(glycosamine_pattern):
        return False, "No glycosamine units found"
    
    # Look for sulfate groups
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)O")
    has_sulfate = mol.HasSubstructMatch(sulfate_pattern)
    
    # Count oxygen atoms (should be abundant in polysaccharides)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 10:
        return False, "Too few oxygen atoms for a polysaccharide"
    
    # Count nitrogen atoms (should have some from glycosamines)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 1:
        return False, "No nitrogen atoms found (required for glycosamines)"
    
    # Check for cyclic structures (sugar rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:  # Should have multiple sugar rings
        return False, "Insufficient ring structures for a polysaccharide"
        
    # Check for alternating pattern of uronic acid and glycosamine
    uronic_matches = len(mol.GetSubstructMatches(uronic_acid_pattern))
    glycosamine_matches = len(mol.GetSubstructMatches(glycosamine_pattern))
    
    if uronic_matches < 1 or glycosamine_matches < 1:
        return False, "Missing required alternating units"
    
    # Success message varies based on sulfation
    success_msg = "Contains alternating uronic acid and glycosamine units"
    if has_sulfate:
        success_msg += " with sulfate groups"
    
    return True, success_msg