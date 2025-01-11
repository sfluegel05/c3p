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

    # Indole pattern - flexible to allow for substitutions that might not disrupt basic indole skeleton
    indole_pattern = Chem.MolFromSmarts("c1[nH]c2ccccc2c1")  # A simpler but broader indole skeleton pattern
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole skeleton found"

    # Count all nitrogen atoms for the alkaloid classification
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    nitrogen_count = len(nitrogen_atoms)
    if nitrogen_count < 2:  # Indole has 1 nitrogen, thus need at least 1 more for alkaloid
        return False, "Not enough nitrogen atoms to be an alkaloid"

    # Check for basic nitrogen groups typical in alkaloids; we assume a structure with at least one more basic nitrogen
    basic_nitrogen_count = sum(1 for atom in nitrogen_atoms if atom.GetFormalCharge() == 0 and atom.GetHybridization() not in [Chem.rdchem.HybridizationType.SP2, Chem.rdchem.HybridizationType.SP])
    if basic_nitrogen_count < 1:
        return False, "No basic nitrogen atom found, atypical for an alkaloid"

    return True, "Contains indole skeleton and typical nitrogen features of alkaloids"