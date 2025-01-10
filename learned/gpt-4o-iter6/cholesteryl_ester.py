"""
Classifies: CHEBI:17002 cholesteryl ester
"""
from rdkit import Chem

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    A cholesteryl ester is defined by the condensation of the carboxy group of any carboxylic acid
    with the 3-hydroxy group of cholesterol.

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

    # Cholesterol structure SMARTS pattern
    cholesterol_pattern = Chem.MolFromSmarts("[C@H]1(C)CC[C@]2(C)CC[C@H]3[C@@H]4CCC=C5C[C@@H](O)CC[C@@]5(C)[C@]4(C)CC[C@]3(C)[C@H]2C1")
    if not mol.HasSubstructMatch(cholesterol_pattern):
        return False, "No cholesterol backbone found"

    # Ester linkage pattern: searching for -O-C(=O)-R where R is a long chain
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][C@H]1C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester linkage found to cholesterol"
    
    return True, "Contains cholesterol backbone with ester linkage, indicating a cholesteryl ester"

# Test the function with a sample SMILES for cholesteryl ester
smiles_example = "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@H](CC[C@]4(C)[C@H]3CC[C@]12C)OC(=O)CCCCCCC\C=C/C\C=C/CCCCC"
result = is_cholesteryl_ester(smiles_example)
print(result)