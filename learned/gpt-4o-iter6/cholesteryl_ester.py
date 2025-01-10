"""
Classifies: CHEBI:17002 cholesteryl ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Broad cholesterol pattern focusing on steroid structure (four rings)
    cholesterol_pattern = Chem.MolFromSmarts("C1CC2C3C(C1)CC(C2)CC4(C3=CC=C4)C")
    if not mol.HasSubstructMatch(cholesterol_pattern):
        return (False, "No cholesterol framework identified")

    # Check for an ester group (-OC(=O)-)
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    if not mol.HasSubstructMatch(ester_pattern):
        return (False, "No ester linkage found")

    # Look for ester linkage at 3β-position in steroids using a less strict pattern
    c3_position_pattern = Chem.MolFromSmarts("OC(=O)[C@H]1CC[C@H]2[C@@H]3CC[C@H]4([C@@H]3CC=C4[C@H]2CC1)C")
    if not mol.HasSubstructMatch(c3_position_pattern):
        return (False, "Ester linkage not correctly connected to the 3β-position of cholesterol")

    return (True, "SMILES indicates a cholesteryl ester structure")

# Example test
smiles_example = "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@H](CC[C@]4(C)[C@@H]3CC=C2C1)OC(=O)CCCCCCC\C=C/C\C=C/CCCCC"
result = is_cholesteryl_ester(smiles_example)
print(result)