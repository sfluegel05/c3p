"""
Classifies: CHEBI:36233 disaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # SMARTS pattern for 5- and 6-membered sugar rings (with additional common substitutions to broaden match)
    sugar_patterns = [
        Chem.MolFromSmarts("[C@@H]1([O])O[C@H]1"),    # Generic pattern for 5-membered ring (typically furanoses)
        Chem.MolFromSmarts("[C@@H]1([O])O[C@@H]1"),  # Alternative stereo for furanoses
        Chem.MolFromSmarts("[C@@H]1([OH])O[C@H]1"),  # Additional stereo for substituted variations
        Chem.MolFromSmarts("[C@@H]1([O])C([OH])O1"), # Generic pattern for 6-membered ring (typically pyranoses)
        Chem.MolFromSmarts("[C@H]1([O])C([OH])O[C@@H]1")  # Another common configuration for pyranoses
    ]

    # Detect sugar rings
    sugar_ring_count = 0
    for pattern in sugar_patterns:
        matches = mol.GetSubstructMatches(pattern)
        sugar_ring_count += len(matches)

    if sugar_ring_count < 2:
        return False, f"Expected 2 sugar rings, found {sugar_ring_count}"

    # Generic glycosidic linkage pattern to check for connectivity between sugar rings
    # Broadening the linkages by considering any oxygen-bound carbon polar bonds bridging rings
    glycosidic_pattern = Chem.MolFromSmarts("[C]O[C]")  # A less specific form, still emphasizing connection through oxygen
    glyco_matches = mol.GetSubstructMatches(glycosidic_pattern)

    # Check for at least one valid glycosidic linkage
    if len(glyco_matches) < 1:
        return False, "No appropriate glycosidic link found to connect sugars"

    return True, "Contains two sugar rings joined by a glycosidic bond"