"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion based on its SMILES string.
    A quaternary ammonium ion has a nitrogen atom with four organyl groups and a positive charge.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quaternary ammonium ion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the quaternary ammonium ion SMARTS pattern
    quaternary_ammonium_pattern = Chem.MolFromSmarts("[N+](C)(C)(C)C")
    
    # Check if the molecule contains a quaternary ammonium substructure
    if mol.HasSubstructMatch(quaternary_ammonium_pattern):
        return True, "Contains a quaternary ammonium ion structure"

    return False, "No quaternary ammonium ion structure found"

# Example usage
example_smiles = "CCCC[N+](CCCC)(CCCC)CCCC"
is_quat_ion, reason = is_quaternary_ammonium_ion(example_smiles)
print(f"Is Quaternary Ammonium Ion: {is_quat_ion}, Reason: {reason}")