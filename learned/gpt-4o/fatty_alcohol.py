"""
Classifies: CHEBI:24026 fatty alcohol
"""
from rdkit import Chem

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    Fatty alcohols have a long carbon chain (3-27+ carbons) and one or more hydroxyl groups attached.

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

    # Check for hydroxyl groups (alcohol) in the molecule
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4;H2,H1]O")  # Generic pattern for a C-OH group
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl groups (alcohol groups) found"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Ensure the main carbon chain is in the range of fatty alcohols
    if c_count < 3:
        return False, "Too few carbon atoms for fatty alcohol"
    elif c_count > 27:
        # Since the definition includes >27, we generally accept such cases
        return True, "Fatty alcohol with chain longer than 27 carbons"

    return True, f"Fatty alcohol with valid carbon chain length: {c_count} carbons"

# Test examples
example_smiles = "CCCCCCCCCCCCCCC(O)CCC"
result, reason = is_fatty_alcohol(example_smiles)
print(result, reason)