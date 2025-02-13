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
    - Core oligosaccharide structure
    - Specific sugar units including KDO and heptose
    - Characteristic fatty acid chains
    - Phosphate groups in many cases
    
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

    # Molecular weight filter
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:  # LPS components are typically large molecules
        return False, "Molecular weight too low for LPS"

    # Count key atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if o_count < 6:  # Need multiple oxygen atoms for glycosidic bonds and hydroxyl groups
        return False, "Too few oxygen atoms for LPS"
    
    if c_count < 15:  # Need sufficient carbons for sugar and lipid components
        return False, "Too few carbon atoms for LPS"

    # Look for sugar patterns more specifically
    pyranose_pattern = Chem.MolFromSmarts("[CR1][OR1][CR1]([OR0])[CR1]([OR0])[CR1]([OR0])[CR1]")
    furanose_pattern = Chem.MolFromSmarts("[CR1][OR1][CR1]([OR0])[CR1]([OR0])[CR1]")
    sugar_matches = len(mol.GetSubstructMatches(pyranose_pattern)) + len(mol.GetSubstructMatches(furanose_pattern))

    # Look for glycosidic linkages more specifically
    glycosidic_pattern = Chem.MolFromSmarts("[CR1][OR1][CR1]")
    glycosidic_matches = len(mol.GetSubstructMatches(glycosidic_pattern))

    # Look for fatty acid chains with hydroxyl groups
    hydroxy_fatty_acid = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][OX2H1]")
    fatty_acid_matches = len(mol.GetSubstructMatches(hydroxy_fatty_acid))

    # Look for carboxylic acid groups
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    acid_matches = len(mol.GetSubstructMatches(acid_pattern))

    # Calculate score based on structural features
    score = 0
    
    # Sugar content
    if sugar_matches >= 2:
        score += 2
    elif sugar_matches == 1:
        score += 1

    # Glycosidic linkages
    if glycosidic_matches >= 2:
        score += 2
    elif glycosidic_matches == 1:
        score += 1

    # Fatty acid components
    if fatty_acid_matches >= 1:
        score += 2

    # Acid groups
    if acid_matches >= 1:
        score += 1

    # Phosphate groups
    if p_count >= 1:
        score += 1

    # Ring count
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    if ring_count >= 3:
        score += 2
    elif ring_count >= 1:
        score += 1

    # Hydroxyl groups (characteristic of sugars)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches >= 4:
        score += 2
    elif hydroxyl_matches >= 2:
        score += 1

    # Final classification
    if score >= 7:  # Require high confidence for classification
        return True, "Contains characteristic lipopolysaccharide features including sugar units and appropriate linkages"
    else:
        return False, f"Insufficient LPS characteristics (score: {score})"