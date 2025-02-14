"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is an oligosaccharide, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for common sugar rings (hexose and pentose sugars, different stereochemistries)
    sugar_patterns = [
        Chem.MolFromSmarts("[C@@H]1(O)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O"), # Pyranose form
        Chem.MolFromSmarts("[C@H]1(O)[C@@H](O)[C@@H](O)O[C@@H]1[C@H](O)"), # Different stereo pattern
        Chem.MolFromSmarts("[C@@H]1([O])[C@H](O)[C@H](O)[C@@H](O)[C@H]1(O)"), # Another stereo variant
    ]
    
    # Define SMARTS pattern for glycosidic linkages
    glycosidic_bond_patterns = [
        Chem.MolFromSmarts("C-O-C1[C@H](O)"), # glycosidic linkage
        Chem.MolFromSmarts("C1OC(O)"), # another example pattern
    ]
    
    # Check for at least two monosaccharide units
    sugar_count = 0
    for pattern in sugar_patterns:
        sugar_count += len(mol.GetSubstructMatches(pattern))
    
    if sugar_count < 2: # Oligosaccharides must have at least two sugar units
        return False, "Insufficient monosaccharide units"
    
    # Check for glycosidic linkages
    glycosidic_count = 0
    for pattern in glycosidic_bond_patterns:
        if mol.HasSubstructMatch(pattern):
            glycosidic_count += 1
    
    if glycosidic_count < 1:
        return False, "No glycosidic linkages found"
    
    return True, "Contains sufficient monosaccharide units and glycosidic linkages indicative of an oligosaccharide"

# This revised version broadens the detection of sugar units and glycosidic linkages to improve identification accuracy for oligosaccharides.