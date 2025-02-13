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

    # Create SMARTS pattern for the steroid structure
    cholesterol_pattern = Chem.MolFromSmarts("[#6]1[#6]2[#6]3[#6]4[C@H]([#6]C=CC4C3CC(C2C1)OC(=O)C)CCCCC")

    # Check if the molecule has a cholesteryl backbone
    if not mol.HasSubstructMatch(cholesterol_pattern):
        return False, "No cholesteryl backbone found"

    # Create SMARTS pattern for ester linkage (R-COO-R)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    
    # Check for presence of ester linkage in the compound
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"

    return True, "Contains cholesteryl backbone with ester linkage"