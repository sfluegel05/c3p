"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: CHEBI:16852 lipopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    Lipopolysaccharides contain:
    - Trisaccharide repeating unit (two heptose units and octulosonic acid)
    - Oligosaccharide side chains
    - 3-hydroxytetradecanoic acid units
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count key atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if o_count < 4:  # Need multiple oxygen atoms for glycosidic bonds
        return False, "Too few oxygen atoms for lipopolysaccharide"
    
    if c_count < 10:  # Need carbon chains for sugar components
        return False, "Too few carbon atoms for lipopolysaccharide"

    # Look for sugar ring patterns characteristic of LPS
    pyranose_pattern = Chem.MolFromSmarts("[CR1][OR1][CR1][CR1][CR1][CR1]")  # 6-membered sugar ring
    heptose_pattern = Chem.MolFromSmarts("[CR1][OR1][CR1][CR1][CR1][CR1][CR1]")  # 7-membered heptose ring
    sugar_matches = len(mol.GetSubstructMatches(pyranose_pattern))
    heptose_matches = len(mol.GetSubstructMatches(heptose_pattern))
    
    # Look for keto-deoxy-octulosonic acid (KDO) pattern
    kdo_pattern = Chem.MolFromSmarts("[CX4][CX3](=O)[CX4][CX4][CX4][CX4][CX4][CX4]")
    kdo_matches = len(mol.GetSubstructMatches(kdo_pattern))

    # Check for minimum sugar components
    if sugar_matches < 1 and heptose_matches < 1:
        return False, "No sugar rings found"

    # Look for hydroxyl groups (characteristic of sugars)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    if hydroxyl_matches < 2:
        return False, "Too few hydroxyl groups for lipopolysaccharide"

    # Look for 3-hydroxytetradecanoic acid pattern
    hydroxy_fatty_acid = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]")
    fatty_acid_matches = len(mol.GetSubstructMatches(hydroxy_fatty_acid))

    # Look for glycosidic bonds
    glycosidic_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4]")
    glycosidic_matches = len(mol.GetSubstructMatches(glycosidic_pattern))
    
    if glycosidic_matches < 1:
        return False, "No glycosidic linkages found"

    # Check for characteristic LPS features
    lps_features = 0
    if heptose_matches >= 1:
        lps_features += 1
    if kdo_matches >= 1:
        lps_features += 1
    if fatty_acid_matches >= 1:
        lps_features += 1
    if glycosidic_matches >= 2:
        lps_features += 1
    
    # Must have at least 2 characteristic features to be classified as LPS
    if lps_features < 2:
        return False, "Insufficient characteristic LPS features"

    # Count rings to verify oligosaccharide structure
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 1:
        return False, "No rings found"

    return True, "Contains characteristic lipopolysaccharide features including sugar units and appropriate linkages"