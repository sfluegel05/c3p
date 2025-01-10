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
    
    # Core flavonoid: 1-benzopyran ring with an attached phenyl ring on position 2
    flavonoid_core_pattern = Chem.MolFromSmarts("c1cc2c(cc1)c(coc2)c3ccccc3")
    if not mol.HasSubstructMatch(flavonoid_core_pattern):
        return False, "No core flavonoid 1-benzopyran structure with aryl substitution found"
    
    # Check for typical flavonoid functional groups (e.g., -OH, =O)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=O")
    
    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)

    if has_hydroxyl or has_carbonyl:
        return True, "Identified flavonoid structure with typical functional groups"

    return True, "Identified basic flavonoid structure (1-benzopyran core with aryl substitution)"