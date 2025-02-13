"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion based on its SMILES string.
    A quaternary ammonium ion has a nitrogen atom bonded to four univalent organic groups,
    giving it a positive charge.

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

    # Define the SMARTS pattern for quaternary ammonium ion
    quaternary_ammonium_pattern = Chem.MolFromSmarts("[N+](C)(C)(C)C")  # This captures nitrogen bonded to 4 carbons with a positive charge
    if mol.HasSubstructMatch(quaternary_ammonium_pattern):
        return True, "Contains quaternary ammonium ion pattern"
    
    return False, "No quaternary ammonium ion pattern found"