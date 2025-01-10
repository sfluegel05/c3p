"""
Classifies: CHEBI:18254 ribonucleoside
"""
from rdkit import Chem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside is any nucleoside where the sugar component is D-ribose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for beta-D-ribofuranose with correct stereochemistry
    # The pattern matches the D-ribose sugar in furanose form with appropriate stereochemistry
    ribose_smarts = "[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O"
    ribose_pattern = Chem.MolFromSmarts(ribose_smarts)
    if ribose_pattern is None:
        return False, "Invalid ribose SMARTS pattern"

    # Find matches for the ribose sugar
    ribose_matches = mol.GetSubstructMatches(ribose_pattern)
    if not ribose_matches:
        return False, "No D-ribose sugar moiety found"

    # Define SMARTS pattern for beta-N-glycosidic bond
    # The anomeric carbon [C@H]1 in ribose is connected to a nitrogen atom in the nucleobase
    glycosidic_bond_smarts = "[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O[C;H1]-[N]"
    glycosidic_bond_pattern = Chem.MolFromSmarts(glycosidic_bond_smarts)
    if glycosidic_bond_pattern is None:
        return False, "Invalid glycosidic bond SMARTS pattern"

    # Check for the glycosidic bond
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No nucleobase attached via beta-N-glycosidic bond found"

    # Optionally, check for the 2'-hydroxyl group to ensure it's not a deoxyribonucleoside
    # 2'-hydroxyl is on the carbon adjacent to the ring oxygen
    hydroxyl_2prime_smarts = "[C@@H](O)[C@H]1O[C@H](O)[C@@H](O)[C@@H]1O"
    hydroxyl_2prime_pattern = Chem.MolFromSmarts(hydroxyl_2prime_smarts)
    if hydroxyl_2prime_pattern is None:
        return False, "Invalid 2'-hydroxyl SMARTS pattern"

    if not mol.HasSubstructMatch(hydroxyl_2prime_pattern):
        return False, "2'-hydroxyl group not found (may be a deoxyribonucleoside)"

    return True, "Molecule is a ribonucleoside with D-ribose sugar and nucleobase attached via beta-N-glycosidic bond"