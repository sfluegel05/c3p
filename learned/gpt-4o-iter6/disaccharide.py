"""
Classifies: CHEBI:36233 disaccharide
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is defined as two monosaccharides joined by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns for 5- and 6-membered sugar rings
    sugar_patterns = [
        Chem.MolFromSmarts("[C@H]1(O)CO[C@H](O)[C@@H]1O"),  # 5-membered sugar ring
        Chem.MolFromSmarts("[C@H]1(O)CO[C@@H]([C@H](O)[C@H]1O)O"),  # 6-membered sugar ring
    ]

    # Count sugar rings in the molecule
    ring_matches = []
    for pattern in sugar_patterns:
        matches = mol.GetSubstructMatches(pattern)
        ring_matches.extend(matches)

    sugar_ring_count = len(set([atom_idx for match in ring_matches for atom_idx in match]))

    if sugar_ring_count < 2:
        return False, f"Expected at least 2 sugar rings, found {sugar_ring_count // 5 if sugar_ring_count >= 5 else 0}."

    # SMARTS pattern for a glycosidic bond (links anomeric carbon to another carbon through oxygen)
    glycosidic_pattern = Chem.MolFromSmarts("[C@H]1(O)[C@H](O)[C@@H](O)[C@@H](O)O1-[C@H2]O")  # Anomeric connectivity

    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No appropriate glycosidic bond found to connect sugars."

    return True, "Contains at least two sugar rings joined by a glycosidic bond."