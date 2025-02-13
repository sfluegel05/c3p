"""
Classifies: CHEBI:24026 fatty alcohol
"""
from rdkit import Chem

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is an aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for hydroxyl group (OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if c_count < 3:
        return False, "Too few carbon atoms for a fatty alcohol"
    elif c_count > 27 and c_count not in [28, 29, 30, 32]:  # Allow specific exceptions greater than 27
        return False, f"Too many carbon atoms ({c_count} carbon atoms)"

    return True, f"Contains hydroxyl group with {c_count} carbon atoms"