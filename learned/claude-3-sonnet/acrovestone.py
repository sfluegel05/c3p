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

    # Look for isoflavone/flavone core patterns (multiple SMARTS to catch variations)
    core_patterns = [
        # Basic isoflavone core
        "O=C1C=COc2ccc(cc12)",
        # Alternative isoflavone pattern
        "O=C1CC(c2ccccc2)Oc2ccccc12",
        # More flexible pattern for chromene core
        "O=C1C=C(Oc2ccccc12)c1ccccc1",
        # Dihydro form
        "O=C1CC(Oc2ccccc12)c1ccccc1"
    ]
    
    has_core = False
    for pattern in core_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_core = True
            break
    
    if not has_core:
        return False, "No isoflavone/flavone core structure found"

    # Check for various sugar moiety patterns
    sugar_patterns = [
        # Basic pyranose pattern
        "[OX2][CH]1[CH][CH][CH]([CH][CH]1)(O)[CH2][OX2]",
        # Alternative sugar pattern
        "[OX2][CH]1[CH][CH][CH]([CH]1)(O)[CH2][OX2]",
        # Glucuronic acid pattern
        "[OX2][CH]1[CH][CH][CH]([CH][CH]1)(O)C(=O)[OH]",
        # More flexible sugar pattern
        "O[CH]1[CH][CH][CH]([CH]1)CO"
    ]
    
    has_sugar = False
    for pattern in sugar_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_sugar = True
            break

    if not has_sugar:
        return False, "No glycosidic moiety found"

    # Count hydroxyl and methoxy groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    methoxy_pattern = Chem.MolFromSmarts("COc")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    methoxy_count = len(mol.GetSubstructMatches(methoxy_pattern))

    # Look for other characteristic groups
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[OH]")
    glucuronic_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    has_sulfate = mol.HasSubstructMatch(sulfate_pattern)
    has_glucuronic = mol.HasSubstructMatch(glucuronic_pattern)

    # Calculate basic properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    aromatic_rings = len(rdMolDescriptors.CalcAromaticRings(mol))
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Check minimum requirements
    if mol_wt < 300:
        return False, "Molecular weight too low for acrovestone compound"
    
    if aromatic_rings < 2:
        return False, "Insufficient aromatic rings"
    
    if o_count < 5:
        return False, "Insufficient oxygen atoms"

    # Build reason string
    features = []
    features.append("Contains isoflavone/flavone core")
    features.append("Contains glycosidic moiety")
    features.append(f"Has {hydroxyl_count} hydroxyl groups")
    if methoxy_count > 0:
        features.append(f"Has {methoxy_count} methoxy groups")
    if has_sulfate:
        features.append("Contains sulfate group")
    if has_glucuronic:
        features.append("Contains glucuronic acid group")
    features.append(f"Contains {aromatic_rings} aromatic rings")
    features.append(f"Molecular weight: {int(mol_wt)}")

    return True, "; ".join(features)