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

    # More comprehensive sugar patterns
    sugar_patterns = [
        # Basic pyranose ring
        Chem.MolFromSmarts("[CR0,CR1]1[CR0,CR1][CR0,CR1][CR0,CR1][CR0,CR1][OR0,OR1]1"),
        # Furanose ring
        Chem.MolFromSmarts("[CR0,CR1]1[CR0,CR1][CR0,CR1][CR0,CR1][OR0,OR1]1"),
        # KDO-like pattern (3-deoxy sugar)
        Chem.MolFromSmarts("[CR0,CR1]1[CR0,CR1](C(=O))[CR0,CR1][CR0,CR1][CR0,CR1][OR0,OR1]1"),
        # Amino sugar
        Chem.MolFromSmarts("[CR0,CR1]1[CR0,CR1]([NX3])[CR0,CR1][CR0,CR1][CR0,CR1][OR0,OR1]1"),
        # Modified sugar with phosphate
        Chem.MolFromSmarts("[CR0,CR1]1[CR0,CR1][CR0,CR1][CR0,CR1][CR0,CR1]([OR0,OR1]1)[OX2]P(=O)([OX2H,OX1-])[OX2H,OX1-]"),
        # Sugar with sulfate
        Chem.MolFromSmarts("[CR0,CR1]1[CR0,CR1][CR0,CR1][CR0,CR1][CR0,CR1]([OR0,OR1]1)[OX2]S(=O)(=O)[OX2H,OX1-]")
    ]
    
    sugar_matches = []
    for pattern in sugar_patterns:
        if pattern is not None:
            matches = mol.GetSubstructMatches(pattern)
            sugar_matches.extend(matches)
    
    if not sugar_matches:
        return False, "No sugar moieties found"

    # Detect fatty acid chains and modifications
    lipid_patterns = [
        # Long alkyl chain
        Chem.MolFromSmarts("[CH2][CH2][CH2][CH2][CH2][CH2]"),
        # Fatty acid
        Chem.MolFromSmarts("[CX3](=O)[OX2H,OX1-,OX2][CH2][CH2][CH2][CH2]"),
        # Acyl chain with hydroxyl
        Chem.MolFromSmarts("[CH2][CH2][CH2][CH]([OH])[CH2][CH2]"),
        # Branched fatty acid
        Chem.MolFromSmarts("[CH2][CH2][CH2][CH]([CH3])[CH2][CH2]")
    ]
    
    lipid_matches = []
    for pattern in lipid_patterns:
        if pattern is not None:
            matches = mol.GetSubstructMatches(pattern)
            lipid_matches.extend(matches)
            
    if not lipid_matches:
        return False, "No fatty acid chains found"

    # Look for sugar-lipid connections
    connection_patterns = [
        # Ester linkage
        Chem.MolFromSmarts("[CR0,CR1]1[CR0,CR1][CR0,CR1][CR0,CR1][CR0,CR1]([OR0,OR1]1)[OX2]C(=O)[CH2,CH1]"),
        # Amide linkage
        Chem.MolFromSmarts("[CR0,CR1]1[CR0,CR1][NX3]C(=O)[CH2][CR0,CR1][CR0,CR1][OR0,OR1]1"),
        # Phosphodiester linkage
        Chem.MolFromSmarts("[CR0,CR1]1[CR0,CR1][CR0,CR1][CR0,CR1][CR0,CR1]([OR0,OR1]1)[OX2]P(=O)([OX2])[OX2][CH2]"),
        # Direct attachment
        Chem.MolFromSmarts("[CR0,CR1]1[CR0,CR1][CR0,CR1][CR0,CR1][CR0,CR1]([OR0,OR1]1)[OX2,NX3][CH2][CH2]")
    ]

    # Check molecule size and composition
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    
    if num_carbons < 12 or num_oxygens < 4 or mol_weight < 300:
        return False, "Molecule too small to be a saccharolipid"

    # Check for characteristic features
    has_connection = False
    for pattern in connection_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            has_connection = True
            break
            
    if not has_connection:
        # For complex structures, check if there are both sugars and lipids in close proximity
        if len(sugar_matches) >= 1 and len(lipid_matches) >= 1:
            has_connection = True

    if not has_connection:
        return False, "No clear connection between sugar and lipid moieties"

    # Build classification reason
    reason = f"Contains {len(sugar_matches)} sugar rings and {len(lipid_matches)} fatty acid chains"
    if mol.HasSubstructMatch(Chem.MolFromSmarts("P(=O)([O-,OH])([O-,OH])")):
        reason += ", with phosphate groups"
    if mol.HasSubstructMatch(Chem.MolFromSmarts("S(=O)(=O)[O-,OH]")):
        reason += ", with sulfate groups"
    
    return True, reason