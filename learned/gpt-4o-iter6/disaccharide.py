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

    # Extended patterns for ring sugars (5- and 6-membered rings with variations)
    sugar_patterns = [
        Chem.MolFromSmarts("C1OC[C@H](O)[C@@H](O)[C@H]1O"),  # Example: D-glucose pyranose
        Chem.MolFromSmarts("C1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)O1"),  # Example: open and closed sugar rings
        Chem.MolFromSmarts("C1OC(CO)C(O)C1O"),  # Example: D-fructose furanose
    ]
    # Identify potential sugar rings in the molecule
    sugar_ring_count = 0
    for pattern in sugar_patterns:
        matches = mol.GetSubstructMatches(pattern)
        sugar_ring_count += len(matches)

    if sugar_ring_count < 2:
        return False, f"Expected 2 sugar rings, found {sugar_ring_count}"

    # Glycosidic linkage pattern (ensure connectivity between sugars)
    glycosidic_pattern = Chem.MolFromSmarts("[C@H:1]O[C@H:2]")  # Specific linkage pattern with stereochemistry
    glyco_matches = mol.GetSubstructMatches(glycosidic_pattern)

    if len(glyco_matches) < 1:
        return False, "No appropriate glycosidic link found to connect sugars"

    return True, "Contains two sugar rings joined by a glycosidic bond"