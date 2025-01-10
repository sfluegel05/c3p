"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
from rdkit import Chem

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    
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

    # Expanded patterns for triterpenoid backbones
    triterpenoid_patterns = [
        Chem.MolFromSmarts("C1CCC2(C)C(C)C3CCC4C(C)C(CCC5=C(CC[C@@H]6C(C)(C)CC(C)(C)O6)C5)C4(C)C3CCC2C1"),  # oleanane
        Chem.MolFromSmarts("C1CCC2(C)C(C)C3CCC4C(C)C(CCC5=C4C=CC=C5)C3(C)C2C1"),  # ursane
        Chem.MolFromSmarts("C1CCC2(C)(C=CC3C2C(C)CCC4(C)C3CC(C)=C(C)C4"),  # lupane
    ]

    # Ensure patterns are valid
    if not all(pattern for pattern in triterpenoid_patterns):
        return None, "One or more triterpenoid backbone patterns could not be parsed"

    if not any(mol.HasSubstructMatch(pattern) for pattern in triterpenoid_patterns):
        return False, "No triterpenoid backbone found"

    # Improved glycosidic linkage presence
    # Checking for ether linkages that are part of the glycosidic bonds
    glycosidic_patterns = [
        Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1"), # glucose-like ring
        Chem.MolFromSmarts("O[C@@H]1[C@H](O[C@H]2O[C@@H](CO)[C@@H](O)[C@@H](O)[C@@H]2O)[C@@H](O)[C@@H](O)[C@H]1"), # complex sugar component
        Chem.MolFromSmarts("O[C@H]([C@@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O)O[C@@H]1O[C@H]2C[C@@H]1[C@H]1"), # diverse glycosidic linkage in complex saponins
    ]

    # Ensure patterns are valid
    if not all(pattern for pattern in glycosidic_patterns):
        return None, "One or more glycosidic patterns could not be parsed"

    if not any(mol.HasSubstructMatch(pattern) for pattern in glycosidic_patterns):
        return False, "No significant glycosidic linkage found with sugar rings"

    return True, "Contains triterpenoid backbone with glycosidic linkage consistent with saponins"