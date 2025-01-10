"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone typically involves a five-membered lactone ring and
    may contain various levels of saturation and substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is likely a tetrahydrofuranone derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for various tetrahydrofuranone-like structures
    patterns = [Chem.MolFromSmarts(smarts) for smarts in [
        "O1C=CC(=O)C1",      # Unsaturated furanone with double bond
        "O1CC(=O)C=C1",      # Another variation with a double bond
        "O1CCOC1=O",         # Saturated with ester group
        "O1CCC(=O)C1",       # Simple saturated cyclic ester
        "O1C(=O)CC=C1",      # Alpha-beta unsaturated lactone
        "O1CC=C(O)C1",       # Unsaturated, with hydroxyl variations
        "O1C(=O)CCC1"        # Various saturated configurations
    ]]

    # Check if the molecule matches any of the defined tetrahydrofuranone patterns
    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a tetrahydrofuranone-like structure"

    return False, "Does not contain a tetrahydrofuranone-like structure"

# Example usage
# Test the function with a sample SMILES string from known examples
test_smiles = "CCCCCCCC1OC(=O)C(=C)C1C(O)=O"
result, reason = is_tetrahydrofuranone(test_smiles)
print(f"Is Tetrahydrofuranone: {result}, Reason: {reason}")