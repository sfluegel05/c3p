"""
Classifies: CHEBI:33916 aldopentose
"""
from rdkit import Chem

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    An aldopentose is a pentose with a potential aldehyde group at one end.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldopentose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 5:
        return False, f"Expected 5 carbon atoms, found {c_count}"
    
    # Look for aldehyde group (-CHO)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H](=[O])")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No aldehyde group found at one end"

    return True, "Contains 5 carbon atoms and an aldehyde group at one end"

# Example usage:
smiles = "O[C@H]1COC(O)[C@H](O)[C@H]1O"  # Example SMILES for testing
result, reason = is_aldopentose(smiles)
print(f"Is aldopentose: {result}, Reason: {reason}")