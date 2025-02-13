"""
Classifies: CHEBI:36233 disaccharide
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is identified by the presence of two sugar moieties
    connected by a glycosidic bond.
    
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
    
    # Define SMARTS patterns for pyranose and furanose rings
    pyranose_pattern = Chem.MolFromSmarts("C1(CO)OCC(CO)O1") # Simplified
    furanose_pattern = Chem.MolFromSmarts("C1(CO)OCCO1")     # Simplified
    
    # Find matches for these patterns
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)
    
    # Ensure exactly two ring units (either pyranoses or furanoses)
    total_rings = len(pyranose_matches) + len(furanose_matches)
    if total_rings != 2:
        return False, f"Expected 2 monosaccharide units, found {total_rings}"
    
    # Define SMARTS pattern for a glycosidic bond type
    glycosidic_pattern = Chem.MolFromSmarts("C-O-C")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    
    # Ensure there is a glycosidic bond connecting the two sugars
    if len(glycosidic_matches) < 1:
        return False, "No glycosidic bond found connecting sugar units"

    return True, "Contains two monosaccharide units joined by a glycosidic bond"