"""
Classifies: CHEBI:38958 indole alkaloid
"""
from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid contains an indole skeleton and typically additional nitrogen atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for indole skeleton
    indole_pattern = Chem.MolFromSmarts("c1ccc2c(c1)[nH]c2")
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole skeleton found"

    # Check for additional nitrogen atoms indicating alkaloid properties
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if nitrogen_count < 2:
        return False, "Not enough nitrogen atoms to be an alkaloid"

    return True, "Contains indole skeleton and additional nitrogen atoms characteristic of alkaloids"