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

    # Isoflavone/chromone core patterns
    core_patterns = [
        # Basic isoflavone core
        "O=C1C(=COc2ccccc12)c3ccccc3",
        # More general chromone core
        "O=C1C=COc2ccccc12",
        # Alternative pattern for isoflavone
        "O=C1Cc2ccccc2Oc1c3ccccc3"
    ]
    
    has_core = False
    for pattern in core_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_core = True
            break
    
    if not has_core:
        return False, "No isoflavone/chromone core structure found"

    # Sugar and glycosidic patterns
    sugar_patterns = [
        # Pyranose ring
        "O1[C@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O)CO",
        # O-glycosidic bond
        "Oc1ccc([C@H]2O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]2O)cc1",
        # Glucuronic acid
        "O1[C@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O)C(=O)O",
        # Basic sugar ring with hydroxyls
        "O1[CH]([CH]O)[CH]([CH]O)[CH]([CH]1O)CO"
    ]
    
    has_sugar = False
    for pattern in sugar_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_sugar = True
            break

    # Count key features
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    methoxy_pattern = Chem.MolFromSmarts("COc")
    o_glycosidic_pattern = Chem.MolFromSmarts("O[CH]1O[CH][CH][CH][CH][CH]1")
    
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    methoxy_count = len(mol.GetSubstructMatches(methoxy_pattern))
    o_glycosidic_count = len(mol.GetSubstructMatches(o_glycosidic_pattern))
    
    # Calculate basic properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    aromatic_rings = len(mol.GetAromaticAtoms()) // 6  # Approximate number of 6-membered aromatic rings
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Scoring system
    score = 0
    score += 3 if has_core else 0
    score += 2 if has_sugar else 0
    score += min(hydroxyl_count, 3)  # Up to 3 points for hydroxyls
    score += min(methoxy_count, 2)   # Up to 2 points for methoxy groups
    score += 1 if aromatic_rings >= 2 else 0
    score += 1 if o_count >= 5 else 0
    score += 2 if o_glycosidic_count > 0 else 0
    score += 1 if 300 <= mol_wt <= 800 else 0

    # Decision making
    if score < 6:
        return False, "Insufficient structural features for acrovestone classification"

    # Build reason string
    features = []
    if has_core:
        features.append("Contains isoflavone/chromone core")
    if has_sugar:
        features.append("Contains glycosidic moiety")
    if hydroxyl_count > 0:
        features.append(f"Has {hydroxyl_count} hydroxyl groups")
    if methoxy_count > 0:
        features.append(f"Has {methoxy_count} methoxy groups")
    if o_glycosidic_count > 0:
        features.append(f"Has {o_glycosidic_count} O-glycosidic bonds")
    features.append(f"Contains {aromatic_rings} aromatic rings")
    features.append(f"Contains {o_count} oxygen atoms")
    features.append(f"Molecular weight: {int(mol_wt)}")

    return True, "; ".join(features)