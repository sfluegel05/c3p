"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    The molecule should have a five-membered lactone ring that includes an oxygen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a tetrahydrofuranone derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define more comprehensive SMARTS patterns for tetrahydrofuranone-like structures
    # Include alpha-beta unsaturated structures, and variants supporting various substitutions
    patterns = [Chem.MolFromSmarts(smarts) for smarts in [
        "C1C(=O)OC=C1",     # Classical tetrahydrofuranone with double-bond adjacencies
        "C1C(=O)OCC1",      # Saturated tetrahydrofuranone
        "C1COC(=O)C=C1",    # Ensuring variations in connectivity
        "C1COC(=O)CC1",     # Including more flexible lactones
        "O1CC=C(=O)C1"      # Another isomer possibility with slight variance
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