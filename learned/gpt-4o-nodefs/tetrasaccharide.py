"""
Classifies: CHEBI:50126 tetrasaccharide
"""
from rdkit import Chem

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide is a carbohydrate composed of four saccharide units linked by glycosidic bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for common saccharide units and linkages
    glucose_pattern = Chem.MolFromSmarts("[C@H]1([O])[C@H](O)[C@@H](O)[C@H](O)[C@H](O1)")
    mannoside_pattern = Chem.MolFromSmarts("[C@H]1([C@H](O)[C@H](O)[C@@H](O)[C@H](CO)[O]1)")
    galactose_pattern = Chem.MolFromSmarts("[C@H]1([O])[C@@H](O)[C@H](O)[C@H](O)[C@H](O1)")
    
    # Common glycosidic linkage pattern: an ether oxygen bridging two saccharides
    glycosidic_linkage_pattern = Chem.MolFromSmarts("[OH]1[C@](O)([CH])C(O)O1")

    # Initialize counts
    saccharide_count = 0
    linkage_count = 0

    # Check for saccharide units
    if mol.HasSubstructMatch(glucose_pattern):
        saccharide_count += len(mol.GetSubstructMatches(glucose_pattern))
    if mol.HasSubstructMatch(mannoside_pattern):
        saccharide_count += len(mol.GetSubstructMatches(mannoside_pattern))
    if mol.HasSubstructMatch(galactose_pattern):
        saccharide_count += len(mol.GetSubstructMatches(galactose_pattern))

    # Check for glycosidic linkages
    if mol.HasSubstructMatch(glycosidic_linkage_pattern):
        linkage_count = len(mol.GetSubstructMatches(glycosidic_linkage_pattern))

    # Assess if it's a tetrasaccharide
    if saccharide_count >= 4 and linkage_count >= 3:
        return True, "Molecule is a tetrasaccharide with appropriate saccharide units and glycosidic linkages"

    # Provide reason for failure if not a tetrasaccharide
    if saccharide_count < 4:
        return False, f"Less than 4 monosaccharide units found, found {saccharide_count}"
    if linkage_count < 3:
        return False, f"Insufficient glycosidic linkages, found {linkage_count}"

    return False, "Unclassified - unexpected structure"