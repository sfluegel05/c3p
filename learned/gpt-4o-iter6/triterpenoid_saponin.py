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

    # Triterpenoid backbone patterns (pentacyclic, dammarane, lupane, etc.)
    triterpenoid_patterns = [
        Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4(C3)CCC5C4(CC5)C"),   # Generalized pentacyclic
        Chem.MolFromSmarts("C1CCC2CC3CCC4C(C)(C)C5CCC(C12)C45"),      # Oleanane
        Chem.MolFromSmarts("C1CCC2CCCC3C2C4OC(C5)(CC3)CCC45"),        # Ursane
        Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2=CC=C4C3=CC=C5[C@H](C4)C5")  # Dammarane
    ]

    # Ensure at least one triterpenoid backbone is present
    backbone_present = any(mol.HasSubstructMatch(pattern) for pattern in triterpenoid_patterns)
    if not backbone_present:
        return False, "No identifiable triterpenoid backbone found"

    # Glycosidic patterns to detect possible sugar moieties
    glycosidic_patterns = [
        Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@@H](O)[C@H](O)[C@H](O1)"),        # Glucopyranoside
        Chem.MolFromSmarts("O[C@H]1[C@@H](CO)[C@H](O)[C@H](O)[C@H](O1)"),       # Different linkage
        Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O1)")        # Alternative stereochemistry
    ]

    # Check for presence and count of glycosidic linkages
    glycosidic_present = any(mol.HasSubstructMatch(pattern) for pattern in glycosidic_patterns)
    if not glycosidic_present:
        return False, "No glycosidic linkages detected in the structure"

    sugar_count = sum(len(mol.GetSubstructMatches(pattern)) for pattern in glycosidic_patterns)
    if sugar_count < 1:
        return False, "Not enough sugar moieties detected"

    return True, "Contains triterpenoid backbone with sufficient glycosidic linkages"