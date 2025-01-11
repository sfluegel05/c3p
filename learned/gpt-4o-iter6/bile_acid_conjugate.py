"""
Classifies: CHEBI:36249 bile acid conjugate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    A bile acid conjugate typically involves a bile acid structure attached to
    hydrophilic/charged groups such as glycine, taurine, or sulfate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of conjugate patterns
    glycine_pattern = Chem.MolFromSmarts("NCC(=O)O")
    taurine_pattern = Chem.MolFromSmarts("C(CNC)S(=O)(=O)O")
    sulfate_pattern = Chem.MolFromSmarts("[O-]S(=O)(=O)[O-]")

    # Match the patterns for conjugates
    has_glycine = mol.HasSubstructMatch(glycine_pattern)
    has_taurine = mol.HasSubstructMatch(taurine_pattern)
    has_sulfate = mol.HasSubstructMatch(sulfate_pattern)

    if not any([has_glycine, has_taurine, has_sulfate]):
        return False, "No known bile acid conjugate pattern found"

    # Check for core bile acid ring structure (cyclopenta[a]phenanthrene)
    bile_acid_core_pattern = Chem.MolFromSmarts("C1CCC2C3[C@@H]4CC[C@]5(C3C[C@@H]2C1)C(=CC=C5)C4")
    if not mol.HasSubstructMatch(bile_acid_core_pattern):
        return False, "No bile acid core structure found"

    # If both the conjugate pattern and the bile acid core are matched, classify as bile acid conjugate
    return True, "Matches bile acid core with conjugate pattern"

# A test run with a bile acid conjugate SMILES
example_smiles = "C[C@H](CCC(=O)NCC(O)=O)[C@H]1CC[C@H]2[C@@H]3[C@H](O)C[C@@H]4C[C@H](O)CC[C@]4(C)[C@H]3C[C@H](O)[C@]12C"
result, reason = is_bile_acid_conjugate(example_smiles)
print(result, reason)