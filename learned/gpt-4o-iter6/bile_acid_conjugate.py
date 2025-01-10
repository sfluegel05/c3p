"""
Classifies: CHEBI:36249 bile acid conjugate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    A bile acid conjugate typically involves a bile acid structure attached to
    hydrophilic/charged groups such as glycine, taurine, sulfate, etc.

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

    # More generalized bile acid core structure
    # Steroid backbone: cyclopenta[a]phenanthrene framework
    steroid_backbone_pattern = Chem.MolFromSmarts("C1(C2CCC3C4C5)[C@@H]5CCC4CC[C@H]3C2C[C@@H]1")  # Updated pattern
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No bile acid core structure found"

    # Known conjugate patterns (expanded or corrected for breadth)
    # Glycine attached through an amide bond
    glycine_pattern = Chem.MolFromSmarts("NCC(=O)O")
    # Taurine attached through an amide bond
    taurine_pattern = Chem.MolFromSmarts("S(=O)(=O)[O,N]CNC(=O)")
    # Sulfate ester
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[O-]")
    # Glucuronic acid
    glucuronic_acid_pattern = Chem.MolFromSmarts("O=C(O[C3H6]*)C(O)C(O)C(O)CO")

    # Match any of the conjugate patterns
    has_glycine = mol.HasSubstructMatch(glycine_pattern)
    has_taurine = mol.HasSubstructMatch(taurine_pattern)
    has_sulfate = mol.HasSubstructMatch(sulfate_pattern)
    has_glucuronic_acid = mol.HasSubstructMatch(glucuronic_acid_pattern)

    if not any([has_glycine, has_taurine, has_sulfate, has_glucuronic_acid]):
        return False, "No known bile acid conjugate pattern found"

    # If both core and conjugate are present
    return True, "Matches bile acid core with conjugate pattern"

# A test with an example SMILES
example_smiles = "C[C@H](CCC(=O)NCC(O)=O)[C@H]1CC[C@H]2[C@@H]3[C@H](O)C[C@@H]4C[C@H](O)CC[C@]4(C)[C@H]3C[C@H](O)[C@]12C"
result, reason = is_bile_acid_conjugate(example_smiles)
print(result, reason)