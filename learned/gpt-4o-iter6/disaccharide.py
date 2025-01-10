"""
Classifies: CHEBI:36233 disaccharide
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is two monosaccharides joined by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for identifying sugar-like structures (generic)
    sugar_patterns = [
        Chem.MolFromSmarts("C1(CO)OC(O)C(O)C1"),  # Generic 5-membered
        Chem.MolFromSmarts("C1(CO)OC(O)C(O)C(O)C1"),  # Generic 6-membered
    ]

    # Detect sugar rings using the SMARTS patterns
    sugar_ring_count = 0
    for pattern in sugar_patterns:
        matches = mol.GetSubstructMatches(pattern)
        sugar_ring_count += len(matches)

    if sugar_ring_count < 2:
        return False, f"Expected 2 sugar rings, found {sugar_ring_count}"

    # Check for any oxygen link between anomeric carbon of one sugar and another (glycosidic bond)
    glycosidic_pattern = Chem.MolFromSmarts("[C@H1]O[C&D2]")  # Basic glycosidic structure
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No appropriate glycosidic bond found to connect sugars"

    return True, "Contains two sugar rings joined by a glycosidic bond"