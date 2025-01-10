"""
Classifies: CHEBI:51963 hopanoid
"""
from rdkit import Chem

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    Hopanoids are a class of pentacyclic triterpenoids with hopane backbones.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern simplifying the hopanoid core structure
    hopanoid_core_pattern = Chem.MolFromSmarts("C1CC[C@]2(C)[C@H](CC[C@@H]3[C@@]4(CC[C@@H]5[C@@]6(C(C)(C)C)CC[C@@H]5C)C)[C@@H](C4)CC3C2)C1")
    
    if hopanoid_core_pattern is None:
        return None, "Failed to create SMARTS pattern"
    
    # Check if the molecule contains the hopanoid core structure
    if mol.HasSubstructMatch(hopanoid_core_pattern):
        return True, "Contains hopanoid-like backbone"

    return False, "No hopanoid backbone found"

# This function attempts a more accurate diagnosis of hopanoids using clarified SMARTS patterns focused on
# the pentacyclic triterpenoid structure, adjusting for detailed stereochemistry and core features.