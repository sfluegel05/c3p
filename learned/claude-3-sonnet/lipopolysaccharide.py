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
    - Core oligosaccharide structure with specific sugars
    - O-antigen repeating units
    - Lipid A component with fatty acid chains
    - Often includes KDO and heptose units
    
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
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    # Basic element requirements
    if o_count < 4:
        return False, "Too few oxygen atoms for LPS"

    # Define characteristic structural patterns
    patterns = {
        # Core oligosaccharide patterns
        'pyranose': '[CR1][OR1][CR1]([OR0])[CR1]([OR0])[CR1]([OR0])[CR1]',
        'kdo_like': '[CR1](=O)[CR1][CR1]([OR0])[CR1]([OR0])[CR1]([OR0])[CR1]',
        'glycosidic': '[CR1][OR1][CR1]',
        'phosphate_ester': '[OR1]P(=O)([OR0])[OR0]',
        'amide': '[NX3][CX3](=O)',
        'hydroxy_acid': '[OX2H1][CX4][CX4][CX3](=O)[OX2]',
        'fatty_chain': '[CX4][CX4][CX4][CX4][CX4][CX4]',
        'acetyl': '[CX3](=O)[CX4][OR0,NX3]'
    }

    # Convert patterns to RDKit molecules
    pattern_mols = {name: Chem.MolFromSmarts(pattern) for name, pattern in patterns.items()}
    
    # Count pattern matches
    matches = {name: len(mol.GetSubstructMatches(pattern)) 
              for name, pattern in pattern_mols.items()}

    # Calculate base score
    score = 0
    
    # Core sugar structure (essential)
    if matches['pyranose'] >= 2:
        score += 3
    elif matches['pyranose'] == 1:
        score += 1

    # KDO-like structure (characteristic)
    if matches['kdo_like'] >= 1:
        score += 2

    # Multiple glycosidic linkages (essential)
    if matches['glycosidic'] >= 3:
        score += 3
    elif matches['glycosidic'] >= 1:
        score += 1

    # Phosphate groups (common in LPS)
    if matches['phosphate_ester'] >= 1:
        score += 2
    
    # Fatty acid components
    if matches['fatty_chain'] >= 2:
        score += 2
    elif matches['fatty_chain'] >= 1:
        score += 1

    # Amide linkages (common in lipid A)
    if matches['amide'] >= 1:
        score += 1

    # Acetyl groups (common modifications)
    if matches['acetyl'] >= 2:
        score += 1

    # Ring count
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count >= 4:
        score += 2
    elif ring_count >= 2:
        score += 1

    # Additional structural requirements
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    # Classification logic
    if score >= 8 and hydroxyl_count >= 3:
        # Must have minimum structural features
        if matches['pyranose'] >= 1 and matches['glycosidic'] >= 2:
            # Check for either phosphate or fatty acid components
            if matches['phosphate_ester'] >= 1 or matches['fatty_chain'] >= 1:
                return True, "Contains characteristic LPS features including core oligosaccharide and lipid components"
    
    # Component classification for smaller molecules
    if 300 <= rdMolDescriptors.CalcExactMolWt(mol) <= 800:
        if (matches['pyranose'] >= 1 and 
            (matches['fatty_chain'] >= 1 or matches['phosphate_ester'] >= 1) and
            hydroxyl_count >= 3):
            return True, "LPS component with characteristic sugar and lipid/phosphate features"
    
    return False, f"Missing essential LPS characteristics (score: {score})"