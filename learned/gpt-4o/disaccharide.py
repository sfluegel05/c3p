"""
Classifies: CHEBI:36233 disaccharide
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide must contain two monosaccharide rings joined by a glycosidic bond.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a disaccharide, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for typical sugar rings with stereochemistry
    pyranose_smarts = "[C@H]1([C@H](O)[C@@H](O)[C@H](O)[C@H]1O)O"
    furanose_smarts = "[C@H]1([C@H](O)[C@H](O)[C@@H]1O)O"

    pyranose_pattern = Chem.MolFromSmarts(pyranose_smarts)
    furanose_pattern = Chem.MolFromSmarts(furanose_smarts)

    # Check for pyranose and furanose rings
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)

    # Expecting exactly two sugar units
    total_sugar_rings = len(pyranose_matches) + len(furanose_matches)
    if total_sugar_rings != 2:
        return False, f"Expected 2 monosaccharide units, found {total_sugar_rings}"

    # Define more specific SMARTS for glycosidic linkage
    glycosidic_smarts = "[CX4H]([OH1,OR0])[OX2H1,OX2H0][CX4H]"
    glycosidic_pattern = Chem.MolFromSmarts(glycosidic_smarts)

    # Check for a glycosidic bond specifically linking the sugar rings
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic bond found connecting sugar units properly"

    return True, "Contains two monosaccharide units joined by a glycosidic bond"