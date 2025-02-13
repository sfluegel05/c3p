"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion based on its SMILES string.

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

    # Define SMARTS pattern for quaternary ammonium ion
    quaternary_ammonium_pattern = Chem.MolFromSmarts("[N+](C)(C)(C)C")
    
    # Check if the molecule has a substructure match for the quaternary ammonium ion
    if mol.HasSubstructMatch(quaternary_ammonium_pattern):
        return True, "Contains a quaternary ammonium ion structure"

    return False, "No quaternary ammonium ion structure found"