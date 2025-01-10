"""
Classifies: CHEBI:25029 leukotriene
"""
"""
Classifies: CHEBI:6051 leukotriene
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    Leukotrienes are icosanoids with conjugated triene system derived from arachidonic acid.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a leukotriene, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H,OX2-,OX2R]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for conjugated triene system characteristic of leukotrienes
    # More specific patterns for conjugated double bonds
    triene_patterns = [
        "[CH2,CH3]-[CH2,CH1]-C=C-C=C-C=C-[CH2,CH1]-[CH2,CH1]",  # Basic pattern
        "[CH2,CH3]-[CH2,CH1]-C=C-C=C-C=C-C(O)-[CH2,CH1]",       # With hydroxyl
        "[CH2,CH3]-[CH2,CH1]-C=C-C=C-C=C-C(S)-[CH2,CH1]"        # With sulfur (cysteinyl)
    ]
    
    has_triene = False
    for pattern in triene_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(patt):
            has_triene = True
            break
            
    if not has_triene:
        return False, "No characteristic conjugated triene system found"
    
    # Check for basic carbon chain with some flexibility
    # Allow for various modifications while maintaining core structure
    chain_pattern = Chem.MolFromSmarts("[CH3,CH2]-[CH2,CH1]-[CH2,CH1]-[CH2,CH1]-[CH2,CH1]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Missing basic carbon chain structure"
    
    # Count carbons to ensure C20 backbone (allowing some variation for derivatives)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18 or c_count > 24:  # Allow some flexibility for derivatives
        return False, "Carbon count outside acceptable range for leukotriene"
    
    # Check for characteristic oxygen pattern
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Insufficient oxygen content for leukotriene"
    
    # Look for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH1]")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    # Check for cysteinyl group (present in cysteinyl leukotrienes)
    cysteinyl_pattern = Chem.MolFromSmarts("SCC(N)C(=O)")
    has_cysteinyl = mol.HasSubstructMatch(cysteinyl_pattern)
    
    # Additional structural check for characteristic spacing between functional groups
    spacing_pattern = Chem.MolFromSmarts("[OH,SH,NH2]-[CH2,CH1]-[CH2,CH1]-[CH2,CH1]-C(=O)[OH]")
    has_spacing = mol.HasSubstructMatch(spacing_pattern)
    
    if not (has_spacing or has_cysteinyl):
        return False, "Missing characteristic functional group spacing"
    
    # Build classification reason
    features = []
    if hydroxyl_count > 0:
        features.append(f"{hydroxyl_count} hydroxyl groups")
    if has_cysteinyl:
        features.append("cysteinyl group")
        
    features_str = ", ".join(features) if features else "basic structure"
    
    return True, f"Leukotriene with conjugated triene system, carboxylic acid, and {features_str}"