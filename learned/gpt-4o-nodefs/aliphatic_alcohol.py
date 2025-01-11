"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    Aliphatic alcohols contain an -OH group attached to a non-aromatic carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for alcohol group pattern (O single bonded to a saturated/sp3 C)
    alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No aliphatic alcohol group found"
        
    # Check for absence of aromatic rings
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Molecule contains aromatic rings"

    return True, "Contains aliphatic alcohol group and no aromatic rings"

# Example usage:
# result, reason = is_aliphatic_alcohol("CCCC(CC)CCC(O)C")
# print(result, reason)