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
    Leukotrienes are icosanoids with conjugated triene system and specific structural features.
    
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
    
    # Check for C20 backbone (allowing some variation for derivatives)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Not a C20 backbone structure"
    
    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for specific leukotriene conjugated triene patterns
    # This pattern looks for three double bonds that could be conjugated
    triene_patterns = [
        Chem.MolFromSmarts("C=CC=CC=C"),  # Basic conjugated triene
        Chem.MolFromSmarts("C=CC=CC=CC"),  # Extended conjugated system
        Chem.MolFromSmarts("CC=CC=CC=CC")  # Alternative conjugated system
    ]
    
    has_triene = False
    for pattern in triene_patterns:
        if mol.HasSubstructMatch(pattern):
            has_triene = True
            break
    
    if not has_triene:
        return False, "No characteristic conjugated triene system found"
    
    # Check for characteristic carbon chain with double bonds
    chain_pattern = Chem.MolFromSmarts("CCCCC[CD2]~[CD2]~[CD2]~[CD2]~[CD2]~[CD2]~[CD2]CCCC")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Missing characteristic carbon chain structure"
    
    # Additional structural features common in leukotrienes
    hydroxyl_pattern = Chem.MolFromSmarts("[OH1]")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    # Check for cysteinyl group (present in cysteinyl leukotrienes)
    cysteinyl_pattern = Chem.MolFromSmarts("SCC(N)C(=O)")
    has_cysteinyl = mol.HasSubstructMatch(cysteinyl_pattern)
    
    # Count rotatable bonds to ensure appropriate flexibility
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Structure too rigid for leukotriene"
    
    # Build classification reason
    features = []
    if hydroxyl_count > 0:
        features.append(f"{hydroxyl_count} hydroxyl groups")
    if has_cysteinyl:
        features.append("cysteinyl group")
    
    # Additional check for characteristic oxygen pattern
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Insufficient oxygen content for leukotriene"
    
    features_str = ", ".join(features) if features else "basic structure"
    
    return True, f"C20 icosanoid with conjugated triene system, carboxylic acid, and {features_str}"