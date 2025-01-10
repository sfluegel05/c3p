"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is a semisynthetic derivative based on its SMILES string.
    Focus on more precise pattern recognition and combination of synthetic features.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule may be a semisynthetic derivative, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns that may suggest semisynthetic modification
    amide_pattern = Chem.MolFromSmarts("C(=O)N")  # Amide group
    ether_pattern = Chem.MolFromSmarts("C-O-C")  # C-O-C pattern
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")  # Carboxylic ester
    complex_ether_pattern = Chem.MolFromSmarts("C-O-C=C")  # Suggestive of specific synthetic ether groups
    chiral_pattern = Chem.MolFromSmarts("[C@]")  # Chiral center indicator
    complex_aromatics_pattern = Chem.MolFromSmarts("c1[cR2][cR2][cR2][cR2][cR2]c1")

    # Count occurrences and combinations
    has_amide = mol.HasSubstructMatch(amide_pattern)
    has_complex_ether = mol.HasSubstructMatch(complex_ether_pattern)
    has_ester = mol.HasSubstructMatch(ester_pattern)
    has_chiral_centers = mol.HasSubstructMatch(chiral_pattern)
    has_complex_aromatics = mol.HasSubstructMatch(complex_aromatics_pattern)

    # Combine indicators: heuristic based threshold for semisynthetic nature
    if sum([has_amide, has_ester, has_complex_ether, has_chiral_centers, has_complex_aromatics]) >= 3:
        return True, "Markers suggesting semisynthetic derivatization are detected"

    # Without explicit synthesis or further context, remains speculative
    return False, "Not enough evidence for semisynthetic derivation"

# Example test cases
print(is_semisynthetic_derivative("CO[C@H]1C[C@H](O[C@H]2[C@H](C)[C@@H](O[C@@H]3O[C@H](C)C[C@@H]([C@H]3OC(C)=O)N(C)C)[C@@H](C)C[C@@]3(CO3)C(=O)[C@H](C)[C@@H](OC(C)=O)[C@@H](C)[C@@H](C)OC(=O)[C@@H]2C)O[C@@H](C)[C@@H]1OC(C)=O"))  # Troleandomycin
print(is_semisynthetic_derivative("CCCCn1c2cc(OC)ccc2c2ccnc(C)c12"))  # N-butylharmine