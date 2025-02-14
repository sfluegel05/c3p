"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
"""
Classifies: CHEBI:35786 aliphatic aldoxime
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    An aliphatic aldoxime is an oxime (C=N-O) derived from an aliphatic aldehyde.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic aldoxime, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for oxime group (C=N-O)
    oxime_pattern = Chem.MolFromSmarts("[CX3]=[NX2][OX2]")
    if not mol.HasSubstructMatch(oxime_pattern):
        return False, "No oxime group found"
    
    # Check for aliphatic carbon chains (at least one aliphatic carbon atom)
    aliphatic_pattern = Chem.MolFromSmarts("[CX4,CX3]")
    if not mol.HasSubstructMatch(aliphatic_pattern):
        return False, "No aliphatic carbon atoms found"
    
    # Check for aromatic rings or cyclic structures (excluding the oxime ring)
    aromatic_pattern = Chem.MolFromSmarts("c")
    ring_pattern = Chem.MolFromSmarts("[R2]")  # Exclude the oxime ring
    if mol.HasSubstructMatch(aromatic_pattern) or mol.HasSubstructMatch(ring_pattern):
        return False, "Contains aromatic rings or non-oxime cyclic structures, not aliphatic"
    
    # Check for aldehyde functional group (-CHO) or aldehyde-like carbon (C=N-O)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]|[CX3]=[NX2][OX2]")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "Not derived from an aliphatic aldehyde"
    
    # Allow for additional oxygen atoms beyond the oxime group
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    return True, "Contains an oxime (C=N-O) group derived from an aliphatic aldehyde"