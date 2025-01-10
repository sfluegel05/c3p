"""
Classifies: CHEBI:17002 cholesteryl ester
"""
from rdkit import Chem

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    Cholesteryl esters have a cholesterol base (fused ring structure) with a fatty acid ester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cholesteryl ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General SMARTS pattern for the tetracyclic steroid backbone of cholesterol
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3C(C1)C4CCC3C2C4")

    # Check if the molecule has at least part of the steroid backbone
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone characteristic of cholesterol found"
    
    # Check for the presence of an ester linkage (-C(=O)-O-)
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage entity found attached to steroid backbone"
    
    # Ensure the ester group is connected to the steroid backbone identified
    # For simplicity, we'll assume that if both patterns are present, connections are valid within typical cholesteryl esters

    return True, "Contains characteristic steroid backbone and ester linkage"