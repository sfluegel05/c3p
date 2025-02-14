"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: phospho sugar
"""
from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is any monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for phosphate ester group attached to an alcohol
    phosphate_ester_pattern = Chem.MolFromSmarts('[OX2H0]-P(=O)([OX1])[OX1]')
    if phosphate_ester_pattern is None:
        return False, "Invalid phosphate ester SMARTS pattern"

    # Search for phosphate ester groups
    phosphate_matches = mol.GetSubstructMatches(phosphate_ester_pattern)
    if not phosphate_matches:
        return False, "No phosphate ester group attached to an alcohol found"

    # Define SMARTS patterns for sugar rings (furanose and pyranose forms)
    sugar_patterns = [
        Chem.MolFromSmarts('C1[C@H]([OH])OC([C@@H]1[OH])[OH]'),  # Furanose ring
        Chem.MolFromSmarts('C1[C@H]([OH])[C@@H]([OH])O[C@@H]([C@H]1[OH])[OH]')  # Pyranose ring
    ]

    sugar_found = False
    for sugar_pattern in sugar_patterns:
        if sugar_pattern is None:
            continue
        if mol.HasSubstructMatch(sugar_pattern):
            sugar_found = True
            break

    if not sugar_found:
        return False, "No monosaccharide moiety found"

    # Check if the phosphate ester is connected to the sugar
    # Find atoms involved in phosphate ester
    phosphate_atoms = set()
    for match in phosphate_matches:
        phosphate_atoms.update(match)

    # Find atoms in the sugar ring
    for sugar_pattern in sugar_patterns:
        sugar_matches = mol.GetSubstructMatches(sugar_pattern)
        for match in sugar_matches:
            sugar_atoms = set(match)
            # Check if any atom is shared between phosphate ester and sugar
            if phosphate_atoms.intersection(sugar_atoms):
                return True, "Contains monosaccharide with an alcoholic hydroxy group esterified with phosphoric acid"

    return False, "Phosphate ester group is not attached to the sugar moiety"