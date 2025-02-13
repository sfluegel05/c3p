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
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for hexose and pentose sugar rings
    hexose_pattern = Chem.MolFromSmarts("[C@H]1(O)[C@@H](O)[C@@H](O)[C@H](O)[C@H]1O")
    pentose_pattern = Chem.MolFromSmarts("[C@@H]1O[C@H](O)[C@H](O)[C@H]1")
    
    # Find sugar rings
    hexose_matches = mol.GetSubstructMatches(hexose_pattern)
    pentose_matches = mol.GetSubstructMatches(pentose_pattern)
    sugar_count = len(hexose_matches) + len(pentose_matches)
    
    if sugar_count < 2:
        return False, f"Found {sugar_count} sugar units, need at least 2 to be an oligosaccharide"

    # Collect sugar atoms
    sugar_atoms = set()
    for match in hexose_matches + pentose_matches:
        sugar_atoms.update(match)
    
    # Define SMARTS pattern for glycosidic linkage (simplified: C-O-C)
    linkage_pattern = Chem.MolFromSmarts("[#6]-O-[#6]")
    linkage_matches = mol.GetSubstructMatches(linkage_pattern)

    # Check if linkages connect different sugar rings
    connected_sugars = any(
        {sugar1, sugar2}.issubset(link)
        for link in linkage_matches for sugar1 in sugar_atoms for sugar2 in sugar_atoms
        if sugar1 != sugar2
    )

    if connected_sugars:
        return True, "Contains sugar units linked by glycosidic bonds"

    return False, "No glycosidic linkages found among sugar units"