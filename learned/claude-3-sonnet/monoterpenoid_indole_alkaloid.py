"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: monoterpenoid indole alkaloids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    These compounds are derived from tryptophan and typically secologanin.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 240 or mol_wt > 1000:
        return False, f"Molecular weight {mol_wt:.1f} outside typical MIA range"

    # Count basic statistics
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    ring_count = rdMolDescriptors.CalcNumRings(mol)

    # Basic requirements
    if n_count == 0:
        return False, "No nitrogen atoms found - required for alkaloid"
    if c_count < 15:
        return False, f"Too few carbons ({c_count}) for MIA structure"
    if ring_count < 3:
        return False, f"Too few rings ({ring_count}) for MIA structure"

    # Look for indole or modified indole cores
    indole_patterns = [
        "c1ccc2[nH]ccc2c1",  # Basic indole
        "c1ccc2nccc2c1",     # Modified indole
        "c1ccc2N=CCc2c1",    # Another modified form
        "c1ccc2NCCc2c1",     # Dihydroindole
        "C1=CC=C2C(=C1)NC=C2",  # Alternative representation
        "C1=CC=C2C(=C1)N=CC2",  # Another variant
        "[#6]1:[#6]:[#6]:[#6]2:[#6]:[#6]:1:[#7]:[#6]:[#6]:2"  # Generic form
    ]
    
    has_indole = False
    for pattern in indole_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_indole = True
            break
            
    if not has_indole:
        return False, "No indole or modified indole core found"

    # Look for characteristic MIA structural features
    mia_scaffold_patterns = [
        # Common MIA ring fusion patterns
        "[#6]1[#6]2[#7][#6][#6]1[#6][#6]2",  # Basic aspidosperma type
        "[#6]1[#6]2[#7][#6][#6]1[#6][#6]2[#8,#7]",  # Oxygenated variant
        "[#6]1[#6]2[#7][#6][#6]([#6]1)[#6][#6]2",  # Strychnos type
        "[#6]1[#6]2[#7][#6][#6]3[#6]1[#6][#6]2[#6]3",  # Iboga type
        "[#6]1[#6]2[#7][#6][#6]([#6]1)[#6][#6]2[#6](=O)",  # Carbonyl variant
        # Bridge patterns
        "[#6]1[#6]2[#6][#7][#6]1[#6][#6]2",  # Common bridge
        "[#6]1[#6]2[#6][#7][#6]1[#6][#6]2[#8]",  # Oxygenated bridge
    ]
    
    mia_pattern_matches = 0
    for pattern in mia_scaffold_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            mia_pattern_matches += 1

    # Check for secologanin-derived features
    secologanin_patterns = [
        "[#6]-[#6](=O)-O[#6]",  # Ester group
        "[#6]=[#6]-[#6]",  # Vinyl group
        "[#6]-[#6](O)-[#6]",  # Alcohol
        "[#6]1-[#6]-[#6]-[#8]-[#6]-1",  # Pyran ring
    ]
    
    secologanin_matches = 0
    for pattern in secologanin_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            secologanin_matches += 1

    # Calculate revised complexity score
    complexity_contributors = [
        ring_count * 2,          # Rings
        n_count * 2,             # Nitrogens
        mia_pattern_matches * 5,  # MIA patterns (increased weight)
        secologanin_matches * 3,  # Secologanin features
        (1 if o_count > 0 else 0) * 2  # Presence of oxygen
    ]
    complexity_score = sum(complexity_contributors)
    
    if complexity_score < 12:
        return False, f"Complexity score ({complexity_score}) too low for typical MIA"

    # Additional structural requirements
    if ring_count > 12 and mia_pattern_matches == 0:
        return False, "Complex structure lacks characteristic MIA patterns"

    # Check C/N ratio (typically between 6-25 for MIAs)
    cn_ratio = c_count / n_count
    if cn_ratio < 6 or cn_ratio > 25:
        return False, f"C/N ratio ({cn_ratio:.1f}) outside typical MIA range"

    features = []
    if mia_pattern_matches > 0:
        features.append(f"{mia_pattern_matches} MIA scaffold patterns")
    if secologanin_matches > 0:
        features.append(f"{secologanin_matches} secologanin-derived features")
    
    feature_str = ", ".join(features) if features else "characteristic structural features"
    return True, f"Contains indole core, complex ring system ({ring_count} rings), and {feature_str}"