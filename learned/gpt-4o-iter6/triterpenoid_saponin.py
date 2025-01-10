"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    A triterpenoid saponin is a terpene glycoside where the terpene moiety is a triterpenoid.

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

    # Refined pattern for triterpenoid backbone: more specific ring patterns often found in triterpenoids
    triterpenoid_pattern = Chem.MolFromSmarts("[C;R]1[C;R][C;R][C;R]2[C;R][C;R][C;R]3[C;R][C;R][C;R]4[C;R][C;R][C;R]5[C;R][C;R][C;R]6[C;R][C;R][C;R]")
    if not mol.HasSubstructMatch(triterpenoid_pattern):
        return False, "No identifiable triterpenoid patterns found"

    # Check for glycosidic linkages, expanded for more common sugar linkages
    glycosidic_patterns = [
        Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@H](O)[C@H](O)[C@H]1O"),  # Glucopyranoside
        Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O"), # Alternate stereochemistry
        Chem.MolFromSmarts("O[C@H]1[C@@H](CO)[C@H](O)[C@H](O)[C@H]1O"),  # Common linkages for sugars
    ]

    glycosidic_present = any(mol.HasSubstructMatch(pattern) for pattern in glycosidic_patterns)
    if not glycosidic_present:
        return False, "No glycosidic linkages detected in structure"

    # Ensure presence of sugar moieties by validating against sugar motifs
    num_sugar_matches = sum(len(mol.GetSubstructMatches(pattern)) for pattern in glycosidic_patterns)
    if num_sugar_matches < 1:
        return False, "Insufficient number of sugar moieties detected"

    return True, "Contains triterpenoid backbone with sufficient glycosidic linkages and sugar moieties"