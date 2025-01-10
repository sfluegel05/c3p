"""
Classifies: CHEBI:17002 cholesteryl ester
"""
from rdkit import Chem

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    Cholesteryl esters have a cholesterol base (a specific tetracyclic ring structure) with a fatty acid ester linkage.

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

    # More detailed SMARTS pattern for the entire tetracyclic steroid backbone of cholesterol
    steroid_pattern = Chem.MolFromSmarts("C1[C@@H]2CC[C@H]3[C@@H](C2)CC=C4C[C@@H](OC(=O))[C@@]3(CC[C@]4(C[C@@H]5[C@@H](CC[C@]6(CCC=C7[C@H]6CCC(=C[C@]7(C[C@@H]5C)C)C)C)C)C)C")

    # Check if the molecule has the characteristic cholesterol steroid backbone
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No cholesterol steroid backbone fully recognized"

    # Check for the presence of a ester linkage connected to the 3-beta OH group
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0][C@@H]1CCC2C1C=CC3C2(C)CCC4=C3C=CC5[C@]4(C)CC[C@@H]5C")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage entity found at the cholesterol's 3-beta position"

    # Confirm that the ester linkage is part of a cholesteryl ester by verifying its connection
    # Assume connections are valid if both patterns are present for cholesteryl esters

    return True, "Contains characteristic cholesterol steroid backbone and ester linkage at correct position"