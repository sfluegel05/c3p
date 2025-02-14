"""
Classifies: CHEBI:35915 sterol ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid with an ester linkage formed by condensation of a carboxylic acid
    with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Simplified steroid skeleton pattern (three 6-membered and one 5-membered ring fused)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3(CCC4)C")
    
    # Simplified ester pattern -C(=O)O-
    ester_pattern = Chem.MolFromSmarts("C(=O)O")

    # Check for steroid skeleton
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid skeleton found"

    # Check for ester bond
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"

    # Further refine check by ensuring ester is on the 3-hydroxy position
    # Currently, not checking specific position due to complexity, but note for future refinements

    return True, "Contains a steroid skeleton with an ester linkage"

# Test the function with a sterol ester example
smiles_example = "CCCCCCCCCCCCCCCCC(=O)O[C@H]1CC[C@]2(C)[C@H]3CC[C@]4(C)[C@H](CC[C@H]4[C@@H]3CC=C2C1)[C@H](C)CCCC(C)C"
print(is_sterol_ester(smiles_example))