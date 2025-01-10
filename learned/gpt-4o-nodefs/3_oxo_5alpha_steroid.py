"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
from rdkit import Chem

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more general SMARTS pattern for 3-oxo and 5alpha-steroid framework
    # Allow flexibility in stereochemistry - important for recognizing variants
    # Note: Ensure SMARTS is broader to capture alternative positions of core the structures
    stereocenters_pattern = "[C@H]1CC[C@H]2[C@H]3CC[C@@H]4"  # Ring centers with stereochemistry
    core_steroid_pattern = Chem.MolFromSmarts(f"{stereocenters_pattern}[C@](C(=O))CC4=C2C1")
    
    # Check if any part of the steroid framework is present
    ring_pattern = Chem.MolFromSmarts("C1C2C3C4CCCC4C3C2C1") # Basic steroid ring system - ignore stereochemistry here

    # Check for the 3-oxo group; separate check for versatility
    oxo_group_pattern = Chem.MolFromSmarts("C(=O)")

    if mol.HasSubstructMatch(ring_pattern) and mol.HasSubstructMatch(oxo_group_pattern):
        if mol.HasSubstructMatch(core_steroid_pattern):
            return True, "Contains 3-oxo group in a compatible 5alpha-steroid framework"
        return False, "Contains 3-oxo group but lacks proper 5alpha-steroid stereocenter configuration"
    return False, "Does not contain 3-oxo group or any steroidal feature"

    return None, None