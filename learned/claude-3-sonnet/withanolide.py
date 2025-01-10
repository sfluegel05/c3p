"""
Classifies: CHEBI:74716 withanolide
"""
"""
Classifies: withanolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide is a C28 steroid lactone with a modified side chain forming 
    a lactone ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons - withanolides typically have around 28 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (24 <= c_count <= 35):  # Range to accommodate derivatives
        return False, f"Carbon count {c_count} not typical for withanolide structure"

    # Multiple steroid core patterns to catch different variations
    steroid_patterns = [
        # Basic steroid core with flexible bonds
        "[#6]1~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~[#6]~1",
        # Alternative pattern with more specific ring connections
        "[#6]1~[#6]~[#6]2~[#6](~[#6]~[#6]1)~[#6]~1~[#6]~[#6]~[#6]3~[#6](~[#6]~[#6]~[#6]3)~[#6]~1~[#6]2",
        # More flexible pattern for modified cores
        "[#6]1~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6]~[#6]~[#6]~[#6]3~[#6]2~[#6]~1"
    ]
    
    core_found = False
    for pattern in steroid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            core_found = True
            break
            
    if not core_found:
        return False, "No steroid-like core structure found"

    # Check for lactone rings - multiple patterns to catch variations
    lactone_patterns = [
        # General δ-lactone pattern
        "O=C1OCC[C@@H]1",
        # More specific withanolide lactone patterns
        "O=C1OC(C)=C(C)C1",
        "O=C1OCC(C)C1",
        # Unsaturated lactone patterns
        "O=C1OC=CC1",
        "O=C1OC(=C)CC1",
        # More general patterns
        "O=C1OC[C,C=]C1",
        "O=C1O[CH2,CH][CH2,CH]C1"
    ]
    
    lactone_found = False
    for pattern in lactone_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat is not None and mol.HasSubstructMatch(pat):
            lactone_found = True
            break
    
    if not lactone_found:
        return False, "No characteristic lactone ring found"

    # Check for oxygen count
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if not (2 <= o_count <= 15):
        return False, f"Oxygen count {o_count} not typical for withanolides"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (350 <= mol_wt <= 900):
        return False, f"Molecular weight {mol_wt} outside typical range for withanolides"

    # Check ring count and ring connectivity
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    if ring_count < 4:
        return False, f"Ring count {ring_count} too low for withanolide structure"

    # Look for characteristic oxygen substitutions (hydroxyl, ketone, or ester groups)
    oxygen_patterns = [
        "[OX2H1]", # hydroxyl
        "[CX3](=[OX1])", # ketone
        "[#6]OC(=O)[#6]" # ester
    ]
    
    oxygen_features = 0
    for pattern in oxygen_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat is not None and mol.HasSubstructMatch(pat):
            oxygen_features += 1
    
    if oxygen_features < 1:
        return False, "Missing characteristic oxygen-containing groups"

    # If all checks pass
    return True, "Matches withanolide structure with steroid core, lactone ring, and characteristic substitution pattern"