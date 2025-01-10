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

    # Create more general SMARTS pattern for the steroid structure (fused tetracyclic ring)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)C3CCC4(C3CC(C2)CC4)OC(=O)")

    # Check if the molecule has a steroid backbone with an ester linked
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No cholesteryl backbone found or improperly attached ester linkage"

    return True, "Contains cholesteryl backbone with ester linkage"