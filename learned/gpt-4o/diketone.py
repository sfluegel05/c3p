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

    # SMARTS pattern for a ketone group
    ketone_pattern = Chem.MolFromSmarts("C(=O)C")

    # Find ketone groups matches
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)

    # Process the matches to count unique ketone groups
    # creating a set to store bonds represented by sorted tuples
    unique_ketone_bonds = set()
    for match in ketone_matches:
        # Sort the carbon-carbon bonded atom indices to ensure set uniqueness
        bond = tuple(sorted((match[0], match[2])))  # ketone is C-O, look at the second carbon (match[2])
        unique_ketone_bonds.add(bond)

    # Check if there are exactly two ketone groups
    num_ketones = len(unique_ketone_bonds)
    if num_ketones == 2:
        return True, "Contains exactly two ketone groups, hence a diketone."
    else:
        return False, f"Contains {num_ketones} ketone groups, not a diketone."

# Examples
result, reason = is_diketone("CCCCCC(=O)CC(=O)CCCCCC")  # SMILES example of a diketone
print(result, reason)

result, reason = is_diketone("CCCCCCC=O")  # Non-diketone
print(result, reason)