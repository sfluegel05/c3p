"""
Classifies: CHEBI:18154 polysaccharide
"""
from rdkit import Chem

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide is composed of more than ten monosaccharide residues linked glycosidically.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None, "Invalid SMILES string")

    # Define adjusted patterns for detecting glycosidic linkages and monosaccharides
    # These patterns better represent actual bonding and structure
    glycosidic_patterns = [
        Chem.MolFromSmarts("OC1OC[C@H](O)[C@@H](O)[C@H]1O"),  # Simplified pyran structure
        Chem.MolFromSmarts("O[C@H]1O[C@@H](CO)[C@@H](O)[C@H](O)[C@H]1"), # Another pyran
    ]
    
    # Detect presence of glycosidic linkages
    has_glycosidic_linkage = any(mol.HasSubstructMatch(pattern) for pattern in glycosidic_patterns)
    if not has_glycosidic_linkage:
        return (False, "No glycosidic linkages found")
    
    # Pattern for monosaccharide detection
    saccharide_patterns = [
        Chem.MolFromSmarts("[C@H]1(O)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O"),  # Symmetrical glucose-like
    ]

    matches = set()
    for pattern in saccharide_patterns:
        matches.update(mol.GetSubstructMatches(pattern))
    
    # Checking if more than 10 monosaccharides are present
    if len(matches) <= 10:
        return (False, f"Only found {len(matches)} monosaccharide units, require more than 10")

    return (True, "Valid polysaccharide containing necessary glycosidic linkages and sufficient monosaccharide units")