"""
Classifies: CHEBI:46640 diketone
"""
from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone contains exactly two ketone (C=O) functionalities.

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

    # SMARTS pattern for a ketone group considering varied environments (e.g., adjacent to O, N, etc.)
    # The ketone part can vary a lot, so we allow the neighboring atom on either side to be wildcard
    ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#8,#7,#6]")  # X3 carbon with a C=O

    # Find ketone groups matches
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)

    # Use the matches list to count unique carbon(kept as number for uniqueness)
    unique_ketone_carbons = set(match[1] for match in ketone_matches)

    # Check if there are exactly two unique ketone carbons
    num_ketone_carbons = len(unique_ketone_carbons)
    if num_ketone_carbons == 2:
        return True, f"Contains exactly two ketone groups, hence a diketone."
    else:
        return False, f"Contains {num_ketone_carbons} ketone groups, not a diketone."

# Examples
result, reason = is_diketone("CCCCCC(=O)CC(=O)CCCCCC")  # SMILES example of a diketone
print(result, reason)

result, reason = is_diketone("CCCCCCC=O")  # Non-diketone
print(result, reason)