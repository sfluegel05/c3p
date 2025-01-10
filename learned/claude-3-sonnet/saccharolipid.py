"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: saccharolipids - lipids containing carbohydrate moieties
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    Saccharolipids are lipids that contain carbohydrate moieties with fatty acids 
    directly attached to the sugar units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sugar rings with more flexible patterns
    # Includes pyranose, furanose, and modified sugars
    sugar_patterns = [
        # Basic pyranose/furanose ring
        Chem.MolFromSmarts("[CR0,CR1]1[CR0,CR1][CR0,CR1][CR0,CR1][CR0,CR1][OR0,OR1]1"),
        # Amino sugar pattern (like in lipid A)
        Chem.MolFromSmarts("[CR0,CR1]1[CR0,CR1]([NH1,NH2])[CR0,CR1][CR0,CR1][CR0,CR1][OR0,OR1]1"),
        # KDO-like pattern
        Chem.MolFromSmarts("[CR0,CR1]1[CR0,CR1](C(=O))[CR0,CR1][CR0,CR1][CR0,CR1][OR0,OR1]1")
    ]
    
    has_sugar = False
    for pattern in sugar_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            has_sugar = True
            break
            
    if not has_sugar:
        return False, "No sugar moieties found"

    # Look for fatty acid chains (more flexible pattern)
    fatty_acid_pattern = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2][CH2]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "No fatty acid chains found"

    # Look for connections between sugars and fatty acids
    # Various types of linkages
    linkage_patterns = [
        # Ester linkage to sugar
        Chem.MolFromSmarts("[CR0,CR1]1[CR0,CR1][CR0,CR1][CR0,CR1][CR0,CR1][OR0,OR1]1-[CH2,CH1,CH0]-O-C(=O)-[CH2]"),
        # Amide linkage (as in lipid A)
        Chem.MolFromSmarts("[CR0,CR1]1[CR0,CR1]-[NH]-C(=O)-[CH2][CH2][CR0,CR1][CR0,CR1][OR0,OR1]1"),
        # Direct attachment
        Chem.MolFromSmarts("[CR0,CR1]1[CR0,CR1][CR0,CR1][CR0,CR1][CR0,CR1][OR0,OR1]1-[CH2]-[CH2][CH2][CH2]")
    ]
    
    has_lipid_sugar_connection = False
    for pattern in linkage_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            has_lipid_sugar_connection = True
            break
            
    if not has_lipid_sugar_connection:
        return False, "No connection between sugar and lipid moieties found"

    # Additional characteristic features
    phosphate_pattern = Chem.MolFromSmarts("[P](=O)([O-,OH])([O-,OH])")
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[O-,OH]")
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern) if phosphate_pattern else False
    has_sulfate = mol.HasSubstructMatch(sulfate_pattern) if sulfate_pattern else False

    # Count key elements for size check
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if num_carbons < 12 or num_oxygens < 4:
        return False, "Molecule too small to be a saccharolipid"

    # Build classification reason
    reason = "Contains sugar moiety with "
    if has_phosphate:
        reason += "phosphate groups, "
    if has_sulfate:
        reason += "sulfate groups, "
    reason += f"fatty acid chains ({len(fatty_acid_matches)} found), "
    reason += "and direct sugar-lipid connections"

    return True, reason