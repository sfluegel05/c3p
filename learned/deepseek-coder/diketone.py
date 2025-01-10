"""
Classifies: CHEBI:46640 diketone
"""
"""
Classifies: CHEBI:50699 diketone
"""
from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone is a compound that contains exactly two ketone functionalities.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the ketone pattern (C=O where the carbon is bonded to two other carbons)
    ketone_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    
    # Find all matches for the ketone pattern
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    # Check if there are exactly two ketone groups
    if len(ketone_matches) != 2:
        return False, f"Found {len(ketone_matches)} ketone groups, need exactly 2"

    # Ensure that the ketone groups are not part of other functional groups
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    acid_pattern = Chem.MolFromSmarts("[OX2H][CX3](=[OX1])")
    
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    
    # If any ketone group is part of an ester, amide, or acid, it's not a diketone
    for match in ketone_matches:
        for ester_match in ester_matches:
            if match[0] == ester_match[1]:  # Check if the carbon of the ketone is part of an ester
                return False, "Ketone group is part of an ester"
        for amide_match in amide_matches:
            if match[0] == amide_match[1]:  # Check if the carbon of the ketone is part of an amide
                return False, "Ketone group is part of an amide"
        for acid_match in acid_matches:
            if match[0] == acid_match[1]:  # Check if the carbon of the ketone is part of an acid
                return False, "Ketone group is part of an acid"

    return True, "Contains exactly two ketone functional groups"