"""
Classifies: CHEBI:2440 acrovestone
"""
"""
Classifies: acrovestone compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_acrovestone(smiles: str):
    """
    Determines if a molecule is an acrovestone compound based on its SMILES string.
    Acrovestone compounds are typically isoflavone derivatives with glycosidic attachments.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acrovestone compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic benzopyran core with carbonyl - more flexible patterns
    core_patterns = [
        # Basic chromone/isoflavone core (more general)
        "O=C1C=COc2cccc12",
        # Alternative pattern allowing for different fusion points
        "O=C1CCOc2cccc12",
        # More flexible pattern for chromene core
        "O1c2ccccc2C(=O)CC1",
        # Pattern catching benzopyran with carbonyl
        "O1c2ccccc2C(=O)C=C1"
    ]
    
    has_core = False
    for pattern in core_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_core = True
            break
    
    if not has_core:
        return False, "No chromone/isoflavone core structure found"

    # More flexible sugar patterns
    sugar_patterns = [
        # Basic pyranose
        "O1[CH][CH][CH][CH]([CH]1)CO",
        # Alternative sugar pattern
        "O1[CH][CH][CH][CH]([CH]1)C",
        # Glucuronic acid
        "O1[CH][CH][CH][CH]([CH]1)C(=O)O",
        # Any cyclic ether with multiple OH groups
        "O1[CH]([CH][CH][CH][CH]1)([CH2]O)"
    ]
    
    has_sugar = False
    for pattern in sugar_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_sugar = True
            break

    # Count key features
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    methoxy_pattern = Chem.MolFromSmarts("COc")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    methoxy_count = len(mol.GetSubstructMatches(methoxy_pattern))
    
    # Calculate basic properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    aromatic_rings = len(rdMolDescriptors.CalcAromaticRings(mol))
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Scoring system
    score = 0
    score += 2 if has_core else 0
    score += 2 if has_sugar else 0
    score += min(hydroxyl_count, 3)  # Up to 3 points for hydroxyls
    score += min(methoxy_count, 2)   # Up to 2 points for methoxy groups
    score += 1 if aromatic_rings >= 2 else 0
    score += 1 if o_count >= 5 else 0
    score += 1 if 300 <= mol_wt <= 800 else 0

    # Decision making
    if score < 5:
        return False, "Insufficient structural features for acrovestone classification"

    # Build reason string
    features = []
    if has_core:
        features.append("Contains chromone/isoflavone core")
    if has_sugar:
        features.append("Contains glycosidic moiety")
    if hydroxyl_count > 0:
        features.append(f"Has {hydroxyl_count} hydroxyl groups")
    if methoxy_count > 0:
        features.append(f"Has {methoxy_count} methoxy groups")
    features.append(f"Contains {aromatic_rings} aromatic rings")
    features.append(f"Contains {o_count} oxygen atoms")
    features.append(f"Molecular weight: {int(mol_wt)}")

    return True, "; ".join(features)