"""
Classifies: CHEBI:50523 butenolide
"""
from rdkit import Chem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide is defined as a gamma-lactone that consists of a 2-furanone skeleton
    and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a broader gamma-lactone pattern for 2-furanone which allows for substituents
    # This SMARTS pattern can capture a 5-membered lactone with variations at positions 3 and 4.
    gamma_lactone_pattern = Chem.MolFromSmarts("O=C1C=COC1")  # Gamma-lactone core
    
    # Check if the structure contains the gamma-lactone pattern
    if not mol.HasSubstructMatch(gamma_lactone_pattern):
        return False, "No gamma-lactone (lactone) structure found or wrong substitutions"

    # Validate the detected substructure with examples
    example_substructures = [
        "O=C1OC=C(C1)",
        "O=C1O/C=C/C1"
    ]
    
    for example in example_substructures:
        example_pattern = Chem.MolFromSmarts(example)
        if mol.HasSubstructMatch(example_pattern):
            return True, f"Structure matches butenolide subpattern: {example}"
    
    return False, "Pattern does not fully match any known butenolide subpattern"

# Example usage
smiles = "O=C1OCC=C1"  # Example of butenolide (but-2-en-4-olide)
result, reason = is_butenolide(smiles)
print(result, reason)