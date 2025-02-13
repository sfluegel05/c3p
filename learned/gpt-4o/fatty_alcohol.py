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
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OH]")  # Ensure carbon-bound OH
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No carbon-bound hydroxyl group found"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Fatty alcohols should ideally have more than 2 carbon atoms
    if c_count < 3:
        return False, "Too few carbon atoms for a fatty alcohol"

    # It's appropriate to evaluate longer chains as fatty alcohols without strict limits >27
    # However, ensure these are majorly carbon chains (i.e., not predominantly aromatic).
    # Assuming anything with more than 90% carbon out of all heavy atoms as aliphatic.
    total_heavy_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    if c_count / total_heavy_atoms < 0.9:
        return False, f"Structure is not primarily aliphatic with {c_count} out of {total_heavy_atoms} heavy atoms being carbon"

    return True, f"Contains carbon-bound hydroxyl group with {c_count} carbon atoms"