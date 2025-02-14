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

    # Define more generic SMARTS patterns for pyranose and furanose sugar rings
    pyranose_smarts = "[C@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@@H]1"
    furanose_smarts = "[C@H]1(O)[C@H](O)[C@H](O)[C@@H]1"
    
    pyranose_pattern = Chem.MolFromSmarts(pyranose_smarts)
    furanose_pattern = Chem.MolFromSmarts(furanose_smarts)

    # Check for pyranose and furanose rings
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)

    # Expecting exactly two sugar units
    total_sugar_rings = len(pyranose_matches) + len(furanose_matches)
    if total_sugar_rings != 2:
        return False, f"Expected 2 monosaccharide units, found {total_sugar_rings}"

    # Define generalized SMARTS for glycosidic linkage
    glycosidic_smarts = "[OD2H1,OD2H0][CX4]([OX2H1,OX2H0])[C]"
    glycosidic_pattern = Chem.MolFromSmarts(glycosidic_smarts)

    # Check for a glycosidic bond specifically linking the sugar rings
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if len(glycosidic_matches) == 0:
        return False, "No glycosidic bond found connecting sugar units properly"

    return True, "Contains two monosaccharide units joined by a glycosidic bond"