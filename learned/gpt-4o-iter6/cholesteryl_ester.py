"""
Classifies: CHEBI:17002 cholesteryl ester
"""
from rdkit import Chem

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    A cholesteryl ester is defined by the ester linkage of any carboxylic acid with
    the 3-hydroxy group of cholesterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cholesteryl ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")

    # More refined cholesterol pattern to recognize key structural features
    cholesterol_pattern = Chem.MolFromSmarts("CC(C)CCC[C@H](C)[C@H]1CC[C@H]2[C@@H](C3=CC=CC=C3)[C@@H]3CC[C@]4([C@H]3CC=C4C2)C1")
    if not mol.HasSubstructMatch(cholesterol_pattern):
        return (False, "No cholesterol framework identified")

    # Check for an ester group (-OC(=O)-)
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    if not mol.HasSubstructMatch(ester_pattern):
        return (False, "No ester linkage found")

    # Verify that the ester is connected to the C3 position of the cholesterol framework
    c3_ester_connection_pattern = Chem.MolFromSmarts("C[C@H](OC(=O))C1CC1")
    if not mol.HasSubstructMatch(c3_ester_connection_pattern):
        return (False, "Ester linkage not correctly connected to C3 hydroxyl of cholesterol")

    return (True, "SMILES indicates a cholesteryl ester structure")

# Example test
smiles_example = "C[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@H](CC[C@]4(C)[C@@]3([H])CC[C@]12C)OC(=O)CCCCCCC\C=C/C\C=C/CCCCC)"
result = is_cholesteryl_ester(smiles_example)
print(result)