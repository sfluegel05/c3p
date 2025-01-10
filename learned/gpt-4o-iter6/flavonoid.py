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

    # Core structure of flavonoid: 1-benzopyran-4-one or chromen-4-one
    flavonoid_core_pattern = Chem.MolFromSmarts("c1cc2oc(=O)ccc2cc1-c3cc(*)ccc3")  # Core and flexible substitution point
    if not mol.HasSubstructMatch(flavonoid_core_pattern):
        return False, "No flavonoid core (1-benzopyran structure) found with aryl substitution"
    
    # Ensure that it has a well-formed aryl substitution on position 2
    substitution_check_pattern = Chem.MolFromSmarts("c1cc2oc(=O)ccc2cc1-c3ccccc3")  # specifically check if an aryl group is present
    if not mol.HasSubstructMatch(substitution_check_pattern):
        return False, "No aryl substitution found at position 2"

    return True, "Identified flavonoid structure with 1-benzopyran core and aryl substitution at position 2"