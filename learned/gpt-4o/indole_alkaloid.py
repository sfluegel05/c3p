"""
Classifies: CHEBI:38958 indole alkaloid
"""
from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid contains an indole skeleton, often with additional nitrogen atoms.

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

    # Look for an indole skeleton, accounting for substitutions and typical connectivity
    indole_pattern = Chem.MolFromSmarts("c1cc2cc[nH]c2c1")
    if not mol.HasSubstructMatch(indole_pattern):
        indole_derivative_pattern = Chem.MolFromSmarts("c1cc2c(c1)[nH]c([nH2])c2")  # Example of a more flexible pattern
        if not mol.HasSubstructMatch(indole_derivative_pattern):
            return False, "No indole skeleton found"

    # Check for additional nitrogen atoms that confer alkaloid properties
    # Count all nitrogen atoms, including non-protonated tertiary amines which are common in alkaloids
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    nitrogen_count = len(nitrogen_atoms)
    if nitrogen_count < 2:
        return False, "Not enough nitrogen atoms to be an alkaloid"

    # Check for basic nitrogen groups typical in alkaloids
    basic_nitrogen_count = sum(1 for atom in nitrogen_atoms if atom.GetFormalCharge() == 0)
    if basic_nitrogen_count < 1:
        return False, "No basic nitrogen atom found, atypical for an alkaloid"

    return True, "Contains indole skeleton and typical nitrogen features of alkaloids"