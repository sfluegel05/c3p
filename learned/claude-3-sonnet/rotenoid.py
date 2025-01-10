"""
Classifies: CHEBI:71543 rotenoid
"""
"""
Classifies: rotenoid compounds
Based on tetrahydrochromenochromene skeleton with specific substitution patterns
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    Rotenoids have a cis-fused tetrahydrochromeno[3,4-b]chromene skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a rotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core structure patterns - multiple patterns to catch different variations
    # Pattern 1: Basic chromenochromene core (more flexible)
    core_pattern1 = Chem.MolFromSmarts("O1c2c([#6]3O[#6]4[#6][#6][#6][#6][#6]4[#6]3=O)cccc2")
    
    # Pattern 2: Alternative core with different saturation
    core_pattern2 = Chem.MolFromSmarts("O1[#6]2[#6]([#6]3O[#6]4[#6][#6][#6][#6][#6]4[#6]3=O)[#6][#6][#6][#6][#6]2")
    
    # Pattern 3: More general pattern for variations
    core_pattern3 = Chem.MolFromSmarts("O1[#6]2[#6]([#6]3O[#6][#6]~[#6]~[#6]~[#6][#6]3=O)~[#6]~[#6]~[#6]~[#6]2")

    if any(pattern is None for pattern in [core_pattern1, core_pattern2, core_pattern3]):
        return False, "Invalid SMARTS patterns"

    # Check ring count (should have at least 4 rings)
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    if ring_count < 4:
        return False, "Insufficient number of rings (need at least 4)"

    # Count oxygen atoms (rotenoids typically have 4+ oxygens)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, f"Too few oxygen atoms ({o_count}), need at least 4"

    # Look for carbonyl group
    carbonyl_pattern = Chem.MolFromSmarts("C=O")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Missing required carbonyl group"

    # Check for core structure
    has_core = False
    core_type = ""
    
    if mol.HasSubstructMatch(core_pattern1):
        has_core = True
        core_type = "standard"
    elif mol.HasSubstructMatch(core_pattern2):
        has_core = True
        core_type = "saturated"
    elif mol.HasSubstructMatch(core_pattern3):
        has_core = True
        core_type = "variant"
        
    if not has_core:
        return False, "Missing chromenochromene core structure"

    # Additional structural features
    features = []
    
    # Check for methoxy groups
    methoxy_pattern = Chem.MolFromSmarts("OC")
    methoxy_count = len(mol.GetSubstructMatches(methoxy_pattern))
    if methoxy_count > 0:
        features.append(f"contains {methoxy_count} methoxy groups")

    # Check for dioxole ring
    dioxole_pattern = Chem.MolFromSmarts("OCO")
    if mol.HasSubstructMatch(dioxole_pattern):
        features.append("contains dioxole ring")

    # Check for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_count > 0:
        features.append(f"contains {hydroxyl_count} hydroxyl groups")

    # Check for characteristic substitution patterns
    if o_count >= 6:
        features.append("highly oxygenated")

    feature_str = ", ".join(features)
    if feature_str:
        feature_str = f" with {feature_str}"

    return True, f"Contains {core_type} chromenochromene core skeleton{feature_str}"