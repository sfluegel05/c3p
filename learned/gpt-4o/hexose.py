"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is defined as a six-carbon monosaccharide with an aldehyde at position 1 (aldohexose)
    or a ketone group at position 2 (ketohexose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check if there are exactly six carbon atoms
    if c_count != 6:
        return False, f"Contains {c_count} carbon atoms, but a hexose requires exactly 6"

    # Define patterns for an aldehyde or ketone
    aldehyde_pattern = Chem.MolFromSmarts("[CH1](=O)")  # Aldehyde C=O
    ketone_pattern = Chem.MolFromSmarts("[CH2][C](=O)[CH2]")  # Ketone C=O in open chain

    # Check for an aldehyde at position 1
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Structure matches an aldohexose (aldehyde group present)"
    
    # Check for a ketone at position 2 or part of a furanose or pyranose
    ketone_patterns = [
        Chem.MolFromSmarts("[C](=O)C"),  # Open chain ketone
        Chem.MolFromSmarts("[C!H0](=O)[C]")  # Cyclic context
    ]
    for pattern in ketone_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Structure matches a ketohexose (ketone group present)"

    # Check for cyclic forms indicative of a furanose or pyranose
    cyclic_pattern = Chem.MolFromSmarts("OC1OC(CC1)C")  # General cyclic hexose form
    if mol.HasSubstructMatch(cyclic_pattern):
        return True, "Structure matches a cyclic hexose form"

    return False, "Does not contain hexose-defining functional groups (aldehyde or ketone)"