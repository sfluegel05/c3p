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

    # Identify carbohydrate rings (common 6-membered sugar rings)
    sugar_ring_pattern_6 = Chem.MolFromSmarts("C1[C@H]([C@H]([C@H]([C@H](O1)O)O)O)O")
    sugar_ring_pattern_5 = Chem.MolFromSmarts("C1[C@H]([C@H]([C@H](O1)O)O)O")
    
    # Find matches for sugar rings
    sugar_ring_matches_6 = mol.GetSubstructMatches(sugar_ring_pattern_6)
    sugar_ring_matches_5 = mol.GetSubstructMatches(sugar_ring_pattern_5)
    
    # Total sugar rings should be two
    total_sugar_rings = len(sugar_ring_matches_6) + len(sugar_ring_matches_5)
    if total_sugar_rings < 2:
        return False, f"Expected 2 sugar rings, found {total_sugar_rings}"

    # Identify glycosidic linkages (C-O-C bridge between two sugars)
    glycosidic_pattern = Chem.MolFromSmarts("C-O-C")
    glyco_matches = mol.GetSubstructMatches(glycosidic_pattern)
    
    # Verify at least one glycosidic linkage
    if len(glyco_matches) < 1:
        return False, "No glycosidic linkage found to connect sugars"

    if total_sugar_rings == 2 and len(glyco_matches) >= 1:
        return True, "Contains two sugar rings joined by a glycosidic bond"
    
    return False, "Failed to confirm disaccharide structure"