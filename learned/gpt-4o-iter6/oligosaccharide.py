"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    An oligosaccharide is a compound with monosaccharide units joined by glycosidic linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an oligosaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns to match typical sugar rings (e.g., hexose, pentose)
    sugar_patterns = [
        Chem.MolFromSmarts("OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1"),  # Hexose
        Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@H](O)[C@H]1"),  # Pentose
        # Add more patterns here for specific sugar alterations if needed
    ]

    # Identify sugar rings
    sugar_count = 0
    sugar_atoms = set()
    for pattern in sugar_patterns:
        matches = mol.GetSubstructMatches(pattern)
        sugar_count += len(matches)
        for match in matches:
            sugar_atoms.update(match)

    if sugar_count < 2:
        return False, f"Found {sugar_count} sugar units, need at least 2 to be an oligosaccharide"

    # Check for glycosidic linkages (C-O-C) between distinct sugar rings
    linkage_pattern = Chem.MolFromSmarts("[CX4]-O-[CX4]")
    linkage_matches = mol.GetSubstructMatches(linkage_pattern)

    connected_sugars = any(
        all(sugar in link for sugar in sugar_atoms)
        for link in linkage_matches
    )

    if connected_sugars:
        return True, "Contains sugar units linked by glycosidic bonds"

    return False, "No glycosidic linkages found among sugar units"