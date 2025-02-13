"""
Classifies: CHEBI:50523 butenolide
"""
from rdkit import Chem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide is defined as a gamma-lactone that consists of a 2-furanone skeleton and its substituted derivatives.

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
    
    # Improved pattern list for gamma-lactone structure:
    gamma_lactone_patterns = [
        Chem.MolFromSmarts("O=C1OC=CC1"),  # Basic 2-furanone core
        Chem.MolFromSmarts("O=C1[O,c]=C[C,c]C1"),  # Allow substitutions on the lactone ring
        Chem.MolFromSmarts("O=C1OC=C([C,c])C1"),  # Allow substitution at various positions 
    ]
    
    # Check if the SMILES structure has the gamma-lactone pattern
    for pattern in gamma_lactone_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Structure matches gamma-lactone (2-furanone) backbone pattern"
    
    return False, "No gamma-lactone (2-furanone) pattern found"

# Example usage
smiles = "O=C1OCC=C1"  # Example of butenolide (but-2-en-4-olide)
result, reason = is_butenolide(smiles)
print(result, reason)