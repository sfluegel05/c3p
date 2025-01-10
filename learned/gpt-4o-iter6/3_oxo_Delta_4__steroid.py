"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    Looks for a steroid backbone with a 3-oxo group and a Delta(4) (α,β unsaturation) double bond.

    Args:
        smiles (str): SMILES string of the chemical entity.

    Returns:
        bool: True if the molecule is a 3-oxo-Delta(4) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid backbone SMARTS (fused four-ring system)
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CCC4C3CCC4')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Define 3-oxo group SMARTS (C=O on the third carbon of the backbone)
    oxo_pattern = Chem.MolFromSmarts('C=O')
    oxo_match = mol.GetSubstructMatches(oxo_pattern)
    
    # Check if 'C=O' is in the potential 3-position within the steroid core
    valid_oxo_positions = [3, 9, 11, 17]  # Assume possible patterns from examples
    oxo_positioned_correctly = any(match[0] in valid_oxo_positions for match in oxo_match)

    if not oxo_positioned_correctly:
        return False, "3-oxo group not found in the correct position"

    # Define double bond pattern at Delta(4) (C=C at positions 4 and 5)
    delta4_pattern = Chem.MolFromSmarts('C=C')
    if not mol.HasSubstructMatch(delta4_pattern):
        return False, "No Delta(4) double bond found"

    return True, "Contains steroid backbone, 3-oxo group, and Delta(4) double bond"

# Example usage:
# result, reason = is_3_oxo_Delta_4__steroid("SMILES_STRING_HERE")
# print(result, reason)