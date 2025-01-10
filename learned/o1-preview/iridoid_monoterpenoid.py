"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.

    Iridoids usually consist of a cyclopentane ring fused to a six-membered oxygen heterocycle
    with specific stereochemistry. Secoiridoids are formed by cleavage of a bond in the cyclopentane ring
    but retain the six-membered oxygen heterocycle and characteristic functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an iridoid monoterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for iridoids and secoiridoids
    # Core iridoid skeleton: cyclopentane fused to a dihydropyran ring (oxygen-containing six-membered ring)
    iridoid_pattern = Chem.MolFromSmarts("""
        [C@H]1[C@H]([#6])[CH2][C@H]2O[C@@H](C=C1)[C@@H]2
        |
        -*:1-:
        """)
    # Secoiridoid pattern: opened cyclopentane ring with preserved pyran ring and characteristic substituents
    secoiridoid_pattern = Chem.MolFromSmarts("""
        [C@H]1[C@H](C=C[O])[C@@H]2O[C@@H](C=C1)[C@@H]2
        |
        -*:1-:
        """)

    # Extended patterns to account for variations
    iridoid_pattern_alt = Chem.MolFromSmarts("""
        [C@H]1[C@@H]2O[C@H](C=C1)[C@H]([C@H](C2)*)*
        |
        -*:1-:
        """)

    # Check for iridoid structure
    if mol.HasSubstructMatch(iridoid_pattern) or mol.HasSubstructMatch(iridoid_pattern_alt):
        return True, "Molecule matches iridoid monoterpenoid core structure"

    # Check for secoiridoid structure
    if mol.HasSubstructMatch(secoiridoid_pattern):
        return True, "Molecule matches secoiridoid monoterpenoid structure"

    # Additional filtering to reduce false positives
    # Check the molecular formula to ensure it's a monoterpenoid (C10 backbone)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    if 'C10' not in formula:
        return False, "Molecule does not have a monoterpenoid carbon skeleton"

    return False, "Molecule does not match iridoid monoterpenoid structural features"