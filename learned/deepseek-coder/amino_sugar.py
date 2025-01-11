"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: CHEBI:37671 amino sugar
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is a sugar with one or more hydroxyl groups replaced by amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ring structure
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No ring structure found"
    
    # Look for at least 2 hydroxyl groups (sugar-like)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Not enough hydroxyl groups for a sugar"

    # Check for amino groups (including N-acetyl and other substitutions)
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(N=O)]")
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    
    # Check for N-acetyl groups specifically
    n_acetyl_pattern = Chem.MolFromSmarts("[NX3;H0]([C&R])(C=O)")
    n_acetyl_matches = mol.GetSubstructMatches(n_acetyl_pattern)
    
    # Count total amino groups (including N-acetyl)
    total_amino = len(amino_matches) + len(n_acetyl_matches)
    if total_amino < 1:
        return False, "No amino groups found"

    # Check if amino groups are attached to ring carbons
    amino_sugar_pattern = Chem.MolFromSmarts("[C&R]-[NX3;H2,H1,H0;!$(N=O)]")
    if not mol.HasSubstructMatch(amino_sugar_pattern):
        return False, "Amino group not attached to ring carbon"

    # More flexible sugar pattern matching
    # Look for at least 3 contiguous carbons with hydroxyl or amino groups
    sugar_pattern = Chem.MolFromSmarts("[C&R][C&R][C&R]")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar-like carbon chain found"

    # Additional check for common sugar functional groups
    sugar_functional_groups = Chem.MolFromSmarts("[C&R][OH]")
    if not mol.HasSubstructMatch(sugar_functional_groups):
        return False, "No sugar-like functional groups found"

    return True, f"Contains a sugar structure with {total_amino} amino group(s) attached"