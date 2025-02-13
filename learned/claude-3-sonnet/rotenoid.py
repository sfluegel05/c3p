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

    # Core structure patterns - more generalized to catch variations
    # Basic chromeno-chromene core (more flexible pattern)
    core_pattern = Chem.MolFromSmarts("O1c2c(C3OC4=CC=CC=C4C3=O)cccc2")
    if core_pattern is None:
        return False, "Invalid core SMARTS pattern"
        
    # Alternative core pattern for different oxidation states
    core_pattern2 = Chem.MolFromSmarts("O1c2c(C3Oc4ccccc4C3=O)cccc2")
    if core_pattern2 is None:
        return False, "Invalid alternative core SMARTS pattern"

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
    if carbonyl_pattern is None:
        return False, "Invalid carbonyl SMARTS pattern"
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Missing required carbonyl group"

    # Check for core structure
    has_core = False
    core_type = ""
    if mol.HasSubstructMatch(core_pattern):
        has_core = True
        core_type = "standard"
    elif mol.HasSubstructMatch(core_pattern2):
        has_core = True
        core_type = "alternative"
        
    if not has_core:
        return False, "Missing chromenochromene core structure"

    # Look for common substituents
    features = []
    
    # Check for methoxy groups (very common in rotenoids)
    methoxy_pattern = Chem.MolFromSmarts("OC")
    if methoxy_pattern is not None:
        methoxy_count = len(mol.GetSubstructMatches(methoxy_pattern))
        if methoxy_count > 0:
            features.append(f"contains {methoxy_count} methoxy groups")

    # Check for dioxole ring (common in some rotenoids)
    dioxole_pattern = Chem.MolFromSmarts("O1COC2=CC=CC12")
    if dioxole_pattern is not None and mol.HasSubstructMatch(dioxole_pattern):
        features.append("contains dioxole ring")

    # Check for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    if hydroxyl_pattern is not None:
        hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
        if hydroxyl_count > 0:
            features.append(f"contains {hydroxyl_count} hydroxyl groups")

    feature_str = ", ".join(features)
    if feature_str:
        feature_str = f" with {feature_str}"

    return True, f"Contains {core_type} chromenochromene core skeleton{feature_str}"