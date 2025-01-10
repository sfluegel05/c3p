"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    A triterpenoid saponin is defined as a terpene glycoside where the terpene moiety is a triterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enhanced pattern for typical triterpenoid scaffolds
    # Representative cyclic patterns: pentacyclic (e.g., oleanane, ursane) and tetracyclic (steroidal)
    pentacyclic_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4(C3)CCC5C4(CC5)C")  # Generalized pentacyclic pattern
    tetracyclic_pattern = Chem.MolFromSmarts("C1CC2C3CC4C(C3)C2CC1C4C")  # Generalized steroidal

    # Check for triterpenoid backbone
    if not (mol.HasSubstructMatch(pentacyclic_pattern) or mol.HasSubstructMatch(tetracyclic_pattern)):
        return False, "No identifiable triterpenoid backbone found"

    # Glycosidic pattern check for glucopyranoside or more complex sugars
    glycosidic_patterns = [
        Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@@H](O)[C@H](O)[C@H](O1)"),  # Glucopyranoside
        Chem.MolFromSmarts("O[C@H]1[C@@H](CO)[C@H](O)[C@H](O)[C@H](O1)"),  # Different linkage
    ]

    glycosidic_present = any(mol.HasSubstructMatch(pattern) for pattern in glycosidic_patterns)
    if not glycosidic_present:
        return False, "No glycosidic linkages detected in the structure"

    # Count detected sugar moieties to ensure glycosidic presence
    sugar_count = sum(len(mol.GetSubstructMatches(pattern)) for pattern in glycosidic_patterns)
    if sugar_count < 1:
        return False, "Not enough sugar moieties detected"

    return True, "Contains triterpenoid backbone with sufficient glycosidic linkages"