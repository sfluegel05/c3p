"""
Classifies: CHEBI:35381 monosaccharide
"""
from rdkit import Chem

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    Monosaccharides are polyhydroxy aldehydes or ketones, with at least three carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")
    
    # Check for aldehyde group (carbon=oxygen with another carbon/hydrogen)
    aldehyde_pattern = Chem.MolFromSmarts("[CH][C]=O")
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)
    
    # Check for ketone group (carbon=oxygen with carbons on both sides)
    ketone_pattern = Chem.MolFromSmarts("[C][C](=[O])[C]")
    has_ketone = mol.HasSubstructMatch(ketone_pattern)
    
    if not (has_aldehyde or has_ketone):
        return (False, "No carbonyl group found (neither aldehyde nor ketone)")
    
    # Count number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return (False, "Too few carbon atoms; monosaccharides must have at least 3")
    
    # Check for the presence of multiple hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[C][O][H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    # Monosaccharides must have multiple hydroxyl groups; count at least 2
    if len(hydroxyl_matches) < 2:
        return (False, f"Insufficient hydroxyl groups found; needed at least 2, found {len(hydroxyl_matches)}")
    
    # Ensure the structure is a single unit (no glycosidic linkages)
    # Monosaccharides should not have ether linkages -[O]- between sugars
    ether_pattern = Chem.MolFromSmarts("[C]-[O]-[C]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if ether_matches:
        return (False, "Structure contains ether linkages, indicating glycosidic linkages")
    
    return (True, "Structure is a monosaccharide with required carbonyl group, hydroxyls, and maintains a single unit structure")