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

    # SMARTS pattern to match typical sugar rings (e.g., hexose, pentose)
    sugar_patterns = [
        Chem.MolFromSmarts("OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1"),  # Hexose
        Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@H](O)[C@H]1"),  # Pentose
        # More patterns per sugar variations (e.g., deoxy sugars) can be added
    ]

    # Look for sugar substructures
    sugar_count = 0
    sugar_atoms = set()
    for pattern in sugar_patterns:
        matches = mol.GetSubstructMatches(pattern)
        sugar_count += len(matches)
        for match in matches:
            sugar_atoms.update(match)

    if sugar_count < 2:
        return False, f"Found {sugar_count} sugar units, need at least 2 to be an oligosaccharide"

    # Check for glycosidic linkages
    linkage_pattern = Chem.MolFromSmarts("[C]-O-[C]")  # Basic C-O-C linkage pattern

    # Ensure linkage connects separate sugar atoms
    for sugar in sugar_atoms:
        env = Chem.FindAtomEnvironmentOfRadiusN(mol, 1, sugar)
        if any(mol.GetBondBetweenAtoms(idx1, idx2).GetSmarts() == "O" for idx1, idx2 in env):
            return True, "Contains sugar units linked by glycosidic bonds"

    return False, "No glycosidic linkages found among sugar units"