"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid is typically a 1-benzopyran structure with an aryl group at position 2,
    but can include various substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Refined pattern for core flavonoid: 1-benzopyran structure with an aryl group at position 2
    flavonoid_core_pattern = Chem.MolFromSmarts("c1cc(-c2ccccc2)c2occ(c2)c1")
    if not mol.HasSubstructMatch(flavonoid_core_pattern):
        return False, "No core flavonoid 1-benzopyran structure with aryl substitution found"
    
    # Check for common functional groups like hydroxyls [OH] or carbonyls (=O)
    # Note: Position is not restricted; we're checking for general presence
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    carbonyl_pattern = Chem.MolFromSmarts("C=O")
    
    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)

    if has_hydroxyl or has_carbonyl:
        return True, f"Identified flavonoid structure with {'hydroxyl' if has_hydroxyl else ''}{' and ' if has_hydroxyl and has_carbonyl else ''}{'carbonyl' if has_carbonyl else ''} functional groups"

    return True, "Identified basic flavonoid structure (1-benzopyran core with aryl substitution)"