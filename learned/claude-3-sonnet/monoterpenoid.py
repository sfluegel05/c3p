"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: CHEBI:35189 monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    Monoterpenoids are derived from monoterpenes (C10 skeleton) but may have some
    carbons removed or rearranged, and often contain oxygen-containing functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get basic molecular properties
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Expanded carbon range (7-15) to include substituted monoterpenoids
    if c_count < 7 or c_count > 15:
        return False, f"Carbon count ({c_count}) outside typical range for monoterpenoids (7-15)"
    
    # Adjusted molecular weight range
    if mol_wt < 100 or mol_wt > 350:
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical range for monoterpenoids"

    # Core structural patterns
    patterns = {
        "isopropyl": "[CH3][CH1]([CH3])[CH1,CH2]",
        "gem_dimethyl": "[CH3][C]([CH3])",
        "cyclohexane": "C1CCCCC1",
        "cyclopentane": "C1CCCC1",
        "bridged_bicycle": "C12CCC1CC2",
        "isoprene_unit": "CC(=C)C",
        "menthane_skeleton": "CC1CCC(C(C)C)CC1",
        "carane_skeleton": "C12CC1C(CC2)",
        "pinane_skeleton": "C12CCC(C1(C)C)C2",
    }
    
    structural_features = []
    
    # Check for core patterns
    for name, smarts in patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            structural_features.append(name.replace("_", " "))

    # Oxygen-containing functional groups
    o_patterns = {
        "alcohol": "[OH1]",
        "ether": "[OR0]",
        "ketone": "[CX3](=[OX1])[#6]",
        "ester": "[#6][CX3](=[OX1])[OX2][#6]",
        "aldehyde": "[CX3H1](=O)[#6]",
        "epoxide": "C1OC1",
    }
    
    o_features = []
    for name, smarts in o_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            o_features.append(name)
    
    if o_features:
        structural_features.append("oxygen-containing groups")

    # Ring analysis
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 0:
        structural_features.append(f"{ring_count} ring(s)")
        
    # Analyze ring systems
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(r) for r in ring_info.AtomRings()]
    
    # Score calculation
    score = 0
    score += len(structural_features)  # Base score from structural features
    score += 1 if 5 <= ring_count <= 3 else 0  # Typical ring count
    score += 1 if any(rs in [5,6] for rs in ring_sizes) else 0  # Common ring sizes
    score += 1 if 7 <= c_count <= 12 else 0  # Ideal carbon range
    score += 1 if o_count > 0 else 0  # Presence of oxygen
    
    # Classification
    if score >= 3 and len(structural_features) >= 2:
        features_str = ", ".join(structural_features)
        return True, f"Monoterpenoid features found: {features_str}"
    else:
        return False, "Insufficient characteristic features for monoterpenoid classification"