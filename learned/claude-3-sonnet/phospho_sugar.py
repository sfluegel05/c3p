"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: CHEBI:37600 phospho sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is a monosaccharide with a phosphate group esterified to a hydroxyl group.

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

    # Look for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("[OX2][P](=[O])([O,OH])[O,OH]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Sugar patterns covering various monosaccharide forms
    sugar_patterns = [
        # Pyranose (6-membered ring)
        "[C,O][C]1[C,O][C,O][C,O][C,O]O1",
        # Furanose (5-membered ring)
        "[C,O][C]1[C,O][C,O][C,O]O1",
        # Linear sugar with multiple hydroxyls
        "[C,O][C]([OH,O])[C]([OH,O])[C]([OH,O])[C]([OH,O])[C,O]",
        # Ketose pattern
        "[C,O][C](=O)[C]([OH,O])[C]([OH,O])[C]([OH,O])[C,O]",
        # Aldose pattern
        "[C,O][C]([OH,O])[C]([OH,O])[C]([OH,O])[C](=O)",
        # Modified sugar patterns
        "[C,O][C]1[C]([NH2,NAc])[C,O][C,O][C,O]O1",
        "[C,O][C]1[C]([OH,O])[C]([OH,O])[C]([C,O])[C,O]O1"
    ]

    has_sugar = False
    for pattern in sugar_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt and mol.HasSubstructMatch(patt):
            has_sugar = True
            break

    if not has_sugar:
        return False, "No monosaccharide structure found"

    # Check for phosphate ester connection to sugar
    sugar_phosphate_patterns = [
        # Primary alcohol phosphate
        "[C]([C]1[O,C][C,O][C,O][C,O][C,O]1)[OX2][P](=[O])([O,OH])[O,OH]",
        # Secondary alcohol phosphate
        "[C]1[O,C][C,O][C,O][C,O][C,O]1[OX2][P](=[O])([O,OH])[O,OH]",
        # Linear sugar phosphate
        "[C]([C]([OH,O])[C]([OH,O])[C]([OH,O]))[OX2][P](=[O])([O,OH])[O,OH]"
    ]

    has_sugar_phosphate = False
    for pattern in sugar_phosphate_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt and mol.HasSubstructMatch(patt):
            has_sugar_phosphate = True
            break

    if not has_sugar_phosphate:
        return False, "No phosphate ester connection to sugar moiety found"

    return True, "Contains monosaccharide with phosphate ester"