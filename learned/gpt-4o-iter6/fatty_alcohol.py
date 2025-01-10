"""
Classifies: CHEBI:24026 fatty alcohol
"""
from rdkit import Chem

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is an aliphatic alcohol with a carbon chain of 3 to >27 atoms, which may be saturated/unsaturated and branched/unbranched.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count total carbon atoms in the molecule
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Check for at least one hydroxyl (-OH) group
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")  # Hydroxyl group pattern
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # Check carbon count range
    if carbon_count < 3:
        return False, f"Too few carbon atoms ({carbon_count}), need at least 3"
    # While unspecified max, more than 27 is in the description, no check needed unless an upper threshold is found relevant

    return True, "Contains a carbon chain with a hydroxyl group(s), characteristic of fatty alcohols"

# The function is designed to identify key features that generally define fatty alcohols as per the defined chemical class.