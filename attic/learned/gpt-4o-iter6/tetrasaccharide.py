"""
Classifies: CHEBI:50126 tetrasaccharide
"""
from rdkit import Chem

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide consists of four monosaccharide units linked by glycosidic bonds.

    Args:
        smiles (str): SMILES representation of the molecule

    Returns:
        bool: True if the molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a broad pattern for cyclic monosaccharides
    sugar_patterns = [
        Chem.MolFromSmarts("C1(C(O)C(O)C(O)C(O)O1)"),  # Generic pyranose
        Chem.MolFromSmarts("C1(C(O)C(O)C(O)O1)"),      # Generic furanose
    ]
    
    # Count detected monosaccharide units
    sugar_matches = sum(len(mol.GetSubstructMatches(pattern)) for pattern in sugar_patterns)
    if sugar_matches < 4:
        return False, f"Only found {sugar_matches} saccharide units, need at least 4 for a tetrasaccharide"
    
    # Detect glycosidic linkages
    glycosidic_link_pattern = Chem.MolFromSmarts("[C]-O-[C]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_link_pattern)
    if len(glycosidic_matches) < 3:
        return False, f"Insufficient glycosidic linkages for a tetrasaccharide, found {len(glycosidic_matches)}"

    # Check for the minimum number of rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, f"Insufficient ring structures for expected tetrasaccharide, found {ring_info.NumRings()} rings"
    
    return True, "Structure matches a tetrasaccharide with at least four monosaccharide units linked by glycosidic bonds"