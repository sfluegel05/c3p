"""
Classifies: CHEBI:50126 tetrasaccharide
"""
from rdkit import Chem

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide comprises four monosaccharide units linked by glycosidic bonds.

    Args:
        smiles (str): SMILES representation of the molecule

    Returns:
        bool: True if the molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more generalized pattern for cyclic monosaccharides, allowing variance in stereochemistry
    # and recognizing common monosaccharide structural motifs
    monosaccharide_patterns = [
        Chem.MolFromSmarts("[C@@H]1(O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)O1"),  # Generic 6-membered sugar (pyranose)
        Chem.MolFromSmarts("[C@@H]1(O)[C@H](O)[C@H](O)[C@@H](O)O1"),  # Generic 5-membered sugar (furanose)
    ]
    
    # Count matches for any sugar pattern; consider we're looking for redundancy in substruct matches
    sugar_matches = sum(len(mol.GetSubstructMatches(pattern)) for pattern in monosaccharide_patterns)
    if sugar_matches != 4:
        return False, f"Only found {sugar_matches} saccharide units, need exactly 4 for a tetrasaccharide"
    
    # Identify glycosidic linkages within the saccharide units
    # Use a generalized ether pattern, including variations in connectivity
    glycosidic_linkage_pattern = Chem.MolFromSmarts("[C]-O-[C]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_linkage_pattern)
    if len(glycosidic_matches) < 3:
        return False, f"Insufficient glycosidic linkages for a tetrasaccharide, found {len(glycosidic_matches)}"

    # Adjusted check for hydroxyls to ensure adequate presence without a strict quota, given the structural variance
    hydroxyl_count = Chem.MolFromSmarts("[OH]")
    num_hydroxyls = len(mol.GetSubstructMatches(hydroxyl_count))
    
    if num_hydroxyls < 8:  # Allow for variance, this is an arbitrary chosen lower bound
        return False, f"Probably insufficient hydroxyl groups for a tetrasaccharide, only found {num_hydroxyls}"

    # Ensure the presence of ring structures typical of cyclic monosaccharides
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, f"Insufficient ring structures for expected tetrasaccharide, found {ring_info.NumRings()} rings"
    
    return True, "Structure matches a tetrasaccharide with four monosaccharide units linked by glycosidic bonds"