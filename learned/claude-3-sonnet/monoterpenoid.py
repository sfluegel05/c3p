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
    
    # Core structural patterns for monoterpenoid skeleton
    core_patterns = {
        "isopropyl": "[CH3][CH1]([CH3])[CH1,CH2]",
        "gem_dimethyl": "[CH3][C]([CH3])",
        "isoprene_unit": "[CH3,CH2][C]([CH3,CH2])=[C,c]",
        "menthane_skeleton": "CC1CCC(C(C)C)CC1",
        "pinane_skeleton": "C12CCC(C1(C)C)C2",
        "carane_skeleton": "C12CC1C(CC2)",
        "thujane_skeleton": "CC(C)C12CC1CC2",
        "fenchane_skeleton": "CC1(C)C2CCC1(C)C2",
        "bornane_skeleton": "CC1(C)C2CCC1(C)CC2",
    }
    
    # Additional patterns for common modifications
    mod_patterns = {
        "cyclohexane": "C1CCCCC1",
        "cyclopentane": "C1CCCC1",
        "bridged_bicycle": "C12CCC1CC2",
        "aromatic_ring": "c1ccccc1",
        "conjugated_system": "C=CC=C",
    }
    
    # Functional group patterns
    func_patterns = {
        "alcohol": "[OH1]",
        "ether": "[OR0]",
        "ketone": "[CX3](=[OX1])[#6]",
        "ester": "[#6][CX3](=[OX1])[OX2][#6]",
        "aldehyde": "[CX3H1](=O)[#6]",
        "epoxide": "C1OC1",
        "carboxylic_acid": "[CX3](=O)[OX2H1]",
        "enol": "[OX2H1][#6X3]=[#6X3]",
        "hemiacetal": "[OX2H1][CX4]([OX2])",
        "amino": "[NX3;H2,H1;!$(NC=O)]",
        "amide": "[NX3][CX3](=[OX1])[#6]",
    }
    
    # Check for core patterns
    core_features = []
    for name, smarts in core_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            core_features.append(name.replace("_", " "))

    # Check for modification patterns
    mod_features = []
    for name, smarts in mod_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            mod_features.append(name.replace("_", " "))

    # Check for functional groups
    func_features = []
    for name, smarts in func_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            func_features.append(name)

    # Ring analysis
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(r) for r in ring_info.AtomRings()]
    
    # Score calculation
    score = 0
    
    # Core features are most important
    score += len(core_features) * 2
    
    # Modification patterns add support
    score += len(mod_features)
    
    # Functional groups are common in monoterpenoids
    if func_features:
        score += min(len(func_features), 3)  # Cap at 3 to avoid over-counting
    
    # Ring count consideration
    if 1 <= ring_count <= 4:
        score += 1
    
    # Common ring sizes (5,6) are typical
    if any(rs in [5,6] for rs in ring_sizes):
        score += 1
    
    # Molecular weight consideration (more flexible range)
    if 100 <= mol_wt <= 500:
        score += 1

    # Combine all features for reporting
    all_features = []
    if core_features:
        all_features.extend(core_features)
    if mod_features:
        all_features.extend(mod_features)
    if func_features:
        all_features.append("oxygen/nitrogen-containing groups")
    if ring_count > 0:
        all_features.append(f"{ring_count} ring(s)")

    # Classification
    # Require either strong core features or combination of features with functional groups
    if (len(core_features) >= 2) or (score >= 4 and len(all_features) >= 2):
        features_str = ", ".join(all_features)
        return True, f"Monoterpenoid features found: {features_str}"
    else:
        return False, "Insufficient characteristic features for monoterpenoid classification"