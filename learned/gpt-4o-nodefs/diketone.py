"""
Classifies: CHEBI:46640 diketone
"""
from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone is defined as having exactly two ketone groups (C=O),
    where each carbonyl carbon is bonded to two atoms which can be carbon or other heteroatoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a diketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for ketone pattern, where carbonyl C may have heteroatom neighbors
    # Broadened pattern: allow [#6,#7,#8] neighbors, common in diketone systems
    ketone_pattern = Chem.MolFromSmarts("[#6,#7,#8]C(=O)[#6,#7,#8]")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    # Count the number of ketone groups
    ketone_count = len(ketone_matches)
    
    if ketone_count == 2:
        return True, "Contains exactly two ketone groups"
    else:
        return False, f"Found {ketone_count} ketone groups, need exactly 2"