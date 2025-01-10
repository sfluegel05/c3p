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

    # Improved pattern for typical triterpenoid scaffolds
    triterpenoid_patterns = [
        Chem.MolFromSmarts("[C;R]12[C;R][C;R][C;R]3[C;R][C;R][C;R]4[C;R][C;R][C;R]5[C;R][C;R][C;D2][C;D3]12"),  # Representative sterane skeleton
    ]

    triterpenoid_found = any(mol.HasSubstructMatch(pattern) for pattern in triterpenoid_patterns)
    if not triterpenoid_found:
        return False, "No identifiable triterpenoid backbone found"

    # Check for common glycosidic bonds - focus on sugar moieties
    glycosidic_patterns = [
        Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@@H](O)[C@H](O)[C@H](O1)"),  # Glucopyranoside
        Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@@H](O1)"), # Different stereochemistry
        Chem.MolFromSmarts("O[C@H]1[C@@H](CO)[C@H](O)[C@H](O)[C@H](O1)")
    ]

    glycosidic_present = any(mol.HasSubstructMatch(pattern) for pattern in glycosidic_patterns)
    if not glycosidic_present:
        return False, "No glycosidic linkages detected in structure"

    # Count detected sugar moieties to ensure enough glycosidic presence
    sugar_count = sum(len(mol.GetSubstructMatches(pattern)) for pattern in glycosidic_patterns)
    if sugar_count < 1:
        return False, "Not enough sugar moieties detected"

    return True, "Contains triterpenoid backbone with sufficient glycosidic linkages and sugar moieties"