"""
Classifies: CHEBI:15693 aldose
"""
from rdkit import Chem

def is_aldose(smiles: str) -> (bool, str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    Aldoses are polyhydroxy aldehydes (H[CH(OH)]nC(=O)H, n >= 2) and their intramolecular hemiacetals.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for an aldehyde group (-C=O)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No aldehyde group found"

    # Look for polyhydroxy structure (multiple -OH groups)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, f"Not enough hydroxyl groups found, got {len(hydroxyl_matches)}"

    # Check for cyclic forms indicating a hemiacetal (specific to sugars)
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings() > 0:
        return False, "No rings found, indicating lack of hemiacetal form"

    # Verify the presence of minimum carbon number for backbone
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 3:  # At least two carbons excluding aldehyde carbon
        return False, f"Insufficient carbon backbone, found {carbon_count}"

    # If all checks are passed, classify as aldose
    return True, "Structure is consistent with an aldose"