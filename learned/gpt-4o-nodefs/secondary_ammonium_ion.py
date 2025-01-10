"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for secondary ammonium ion 
    secondary_ammonium_pattern = Chem.MolFromSmarts("[NH2+][!#1;R0][!#1;R0]")  # Two non-hydrogen, non-ring bonds

    # Check for the secondary ammonium ion pattern
    if mol.HasSubstructMatch(secondary_ammonium_pattern):
        return True, "Contains secondary ammonium ion with [NH2+] and two organic substituents"
    else:
        return False, "No secondary ammonium ion structure found"

# Example usage:
# print(is_secondary_ammonium_ion("CC[NH2+]CC"))