"""
Classifies: CHEBI:71543 rotenoid
"""
"""
Classifies: rotenoid compounds
Based on tetrahydrochromeno[3,4-b]chromene skeleton with specific substitution patterns
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

    # Basic requirements first
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings (need at least 4)"

    # Count oxygen atoms (rotenoids typically have 4+ oxygens)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, f"Too few oxygen atoms ({o_count}), need at least 4"

    # Core structure patterns - multiple simpler patterns that must all match
    # Basic chromene pattern
    chromene = Chem.MolFromSmarts("O1CCc2ccccc12")
    
    # Fused ring system with ketone
    fused_system = Chem.MolFromSmarts("O1c2ccccc2C(=O)C2=C1cccc2")
    
    # Alternative pattern for the core
    core_alt = Chem.MolFromSmarts("O1c2c(C3OC4ccccc4C3=O)cccc2")
    
    # Pattern for characteristic oxygen bridge
    oxygen_bridge = Chem.MolFromSmarts("O1CCc2c1cc1OC(=O)c3c1cc2")

    if not all(pat is not None for pat in [chromene, fused_system, core_alt, oxygen_bridge]):
        return False, "Invalid SMARTS patterns"

    # Check for presence of key structural features
    matches = []
    if mol.HasSubstructMatch(chromene):
        matches.append("chromene")
    if mol.HasSubstructMatch(fused_system):
        matches.append("fused ketone system")
    if mol.HasSubstructMatch(core_alt):
        matches.append("core structure")
    if mol.HasSubstructMatch(oxygen_bridge):
        matches.append("oxygen bridge")

    if len(matches) < 2:  # Need at least two confirming patterns
        return False, "Missing key structural features of rotenoid skeleton"

    # Look for characteristic substitution patterns
    features = []
    
    # Check for methoxy groups
    methoxy_pattern = Chem.MolFromSmarts("OC")
    methoxy_count = len(mol.GetSubstructMatches(methoxy_pattern))
    if methoxy_count > 0:
        features.append(f"{methoxy_count} methoxy groups")

    # Check for dioxole ring (common in rotenoids)
    dioxole_pattern = Chem.MolFromSmarts("OCO")
    if mol.HasSubstructMatch(dioxole_pattern):
        features.append("dioxole ring")

    # Check for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_count > 0:
        features.append(f"{hydroxyl_count} hydroxyl groups")

    # Check for ketone group (essential for most rotenoids)
    ketone_pattern = Chem.MolFromSmarts("C(=O)C")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "Missing required ketone group"

    # Combine structural features into description
    feature_str = ", ".join(features)
    if feature_str:
        feature_str = f" with {feature_str}"

    matched_str = ", ".join(matches)
    return True, f"Contains rotenoid skeleton ({matched_str}){feature_str}"