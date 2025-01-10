"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid is characterized by a 1-benzopyran structure with an aryl group at position 2.

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

    # Core structure of flavonoid: 1-benzopyran with flexible aryl substitution
    flavonoid_core_pattern = Chem.MolFromSmarts("c1cc2occ(=O)cc2c1-c3:c4:c:c:c:c:c:c4c:c3") # Captures variations with flexible extensions
    if not mol.HasSubstructMatch(flavonoid_core_pattern):
        return False, "No flavonoid core (1-benzopyran structure) with appropriate substitution found"
    
    # Ensure that it has a well-formed aryl substitution on position 2
    aryl_substitution_check = Chem.MolFromSmarts("c1cc2occ(=O)cc2c1-c3c:cc:cc:c3")  # specifically check if an aryl group is present
    if not mol.HasSubstructMatch(aryl_substitution_check):
        return False, "No aryl substitution found at position 2"
    
    # Additional checks for common flavonoid modifications (e.g., glycosylations)
    glycoside_pattern = Chem.MolFromSmarts("[OX2H1]C[CX4H2][O]")
    if mol.HasSubstructMatch(glycoside_pattern):
        return True, "Identified flavonoid structure with glycosylation"

    return True, "Identified flavonoid structure with 1-benzopyran core and aryl substitution at position 2"