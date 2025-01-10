"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 6 carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons != 6:
        return False, f"Number of carbon atoms is {num_carbons}, not equal to 6"
    
    # Check for hydroxyl groups (-OH)
    num_hydroxyls = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
    if num_hydroxyls < 3: # Hexoses typically have multiple hydroxyls
        return False, f"Insufficient hydroxyl groups, found {num_hydroxyls}"
    
    # Check for a ring structure
    if not mol.GetRingInfo().NumRings():
        return False, "No ring structures detected, hexoses typically form rings"

    # Check for presence of aldehyde or ketone form
    # Aldose pattern: O=C[CH]
    aldehyde_pattern = Chem.MolFromSmarts("O=CO")
    # Ketose pattern: C(=O)C
    ketone_pattern = Chem.MolFromSmarts("C(=O)C")
    if not mol.HasSubstructMatch(aldehyde_pattern) and not mol.HasSubstructMatch(ketone_pattern):
        return False, "Neither aldehyde nor ketone functional groups observed"

    return True, "Structure matches criteria for hexose: 6 carbons, multiple hydroxyl groups, cyclic form"