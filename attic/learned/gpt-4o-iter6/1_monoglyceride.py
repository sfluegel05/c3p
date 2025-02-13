"""
Classifies: CHEBI:35759 1-monoglyceride
"""
from rdkit import Chem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride is defined as a monoglyceride in which the acyl substituent is located at position 1.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern definitively targeting 1-monoglyceride structure:
    # Glycerol backbone with an ester linkage at the primary position (position 1)
    # Primary defined by `[C@H]` for stereochemistry if defined, -O-C(=O)
    pattern1 = Chem.MolFromSmarts('[C@@H](CO)C(=O)O')  # Enhanced specificity for chiral center
    pattern2 = Chem.MolFromSmarts('[C@H](CO)C(=O)O')   # Another possibility based on chiral variations
    pattern3 = Chem.MolFromSmarts('[C](CO)C(=O)O')     # General form without stereochemistry
    
    if mol.HasSubstructMatch(pattern1) or mol.HasSubstructMatch(pattern2) or mol.HasSubstructMatch(pattern3):
        return True, "Contains a glycerol backbone with acyl linkage at position 1"

    return False, "Does not have a glycerol backbone with acyl linkage at position 1"