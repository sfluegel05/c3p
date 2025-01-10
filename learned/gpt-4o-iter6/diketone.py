"""
Classifies: CHEBI:46640 diketone
"""
from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone is defined as a compound that contains two ketone functionalities.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is considered a diketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define ketone group pattern: A carbonyl carbon single-bonded to carbon on at least one side
    # The SMARTS pattern ensures the central carbon is double-bonded to oxygen and single-bonded to at least one carbon
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[C]")
    if ketone_pattern is None:
        return False, "Invalid ketone SMARTS pattern"

    # Find substructure matches for ketone groups
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)

    # Check for exactly two ketone groups, while ensuring non-overlapping matches
    ketone_carbon_atoms = {match[0] for match in ketone_matches}  # store unique carbon atoms of carbonyl groups
    
    if len(ketone_carbon_atoms) == 2:
        return True, "Contains exactly 2 ketone groups, sufficient for diketone classification"
    
    if len(ketone_carbon_atoms) < 2:
        return False, f"Found {len(ketone_carbon_atoms)} ketone groups, not enough for diketone"
    
    return False, f"Found {len(ketone_carbon_atoms)} ketone groups, exactly 2 needed for diketone"