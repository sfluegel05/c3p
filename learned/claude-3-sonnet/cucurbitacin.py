"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Classifies: cucurbitacin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def validate_smarts(smarts_list):
    """Validate SMARTS patterns and return only valid ones"""
    valid_patterns = []
    for smarts in smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is not None:
            valid_patterns.append(pattern)
    return valid_patterns

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids with specific structural features.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic molecular properties
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 25:  # Cucurbitacins typically have 30-32 carbons
        return False, "Too few carbons for cucurbitacin"
    
    if o_count < 4:  # Typically have multiple oxygen-containing groups
        return False, "Too few oxygens for cucurbitacin"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 450 or mol_wt > 1000:  # Typical range for cucurbitacins
        return False, "Molecular weight outside typical range for cucurbitacin"

    # Define and validate SMARTS patterns
    core_patterns = validate_smarts([
        # Basic tetracyclic core patterns with variations
        "C1CC2CCC3C4CCC(C4)C3C2C1",
        "C1CC2CCC3C(CCC4CC(C4)C3)C2C1",
        "C1CC2C(CCC3C4CCCC34)CC2C1"
    ])
    
    if not core_patterns:
        return False, "Failed to create valid structural patterns"

    # Check for core structure
    has_core = any(mol.HasSubstructMatch(pattern) for pattern in core_patterns)
    if not has_core:
        return False, "Missing characteristic tetracyclic core"

    # Functional group patterns
    functional_groups = validate_smarts([
        "[OH]",  # Hydroxyl
        "C(=O)",  # Ketone
        "CC(C)(O)",  # Typical side chain fragment
        "C=CC(=O)",  # α,β-unsaturated ketone
        "OC(=O)",  # Ester
        "C=C",  # Double bond
    ])

    # Ring analysis
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    six_membered_rings = sum(1 for size in ring_sizes if size == 6)
    five_membered_rings = sum(1 for size in ring_sizes if size == 5)

    # Scoring system
    score = 0
    
    # Core structure checks
    score += 2 if has_core else 0
    score += 1 if (six_membered_rings >= 3 and five_membered_rings >= 1) else 0
    
    # Functional group checks
    for pattern in functional_groups:
        matches = len(mol.GetSubstructMatches(pattern))
        if matches > 0:
            score += 1
            if matches >= 2:
                score += 1

    # Molecular property checks
    score += 1 if 450 <= mol_wt <= 800 else 0
    score += 1 if o_count >= 6 else 0
    score += 1 if 30 <= c_count <= 35 else 0

    # Additional structural features
    if ring_info.NumRings() >= 4:
        score += 2

    # Decision threshold
    if score >= 8:
        return True, "Contains cucurbitacin structural features: tetracyclic core, multiple oxygen-containing groups, and characteristic substitution pattern"
    else:
        return False, f"Insufficient cucurbitacin characteristics (score: {score})"