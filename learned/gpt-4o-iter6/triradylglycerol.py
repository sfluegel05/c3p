"""
Classifies: CHEBI:76579 triradylglycerol
"""
from rdkit import Chem

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the glycerol backbone pattern: C-C-C
    glycerol_pattern = Chem.MolFromSmarts("[CH2]C[CH2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Check for ester groups (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    
    # Check for three substituent groups: acyl, alkyl, or alk-1-enyl
    if ester_matches != 3:
        return False, f"Found {ester_matches} ester groups, need 3 for classification"
    
    # If necessary, additional checks for the substituent type can be added here
    # to distinguish between acyl, alkyl, or alk-1-enyl.

    return True, "Contains glycerol backbone with three ester groups"

# Example usage:
# result, reason = is_triradylglycerol("Your SMILES string here")
# print(result, reason)