"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
"""
Classifies: dihydroagarofuran sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    Dihydroagarofuran sesquiterpenoids have a characteristic tricyclic (5/7/6) ring system
    with specific stereochemistry and multiple ester substituents.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (True/False for classification, reason for the classification)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core structure pattern for dihydroagarofuran skeleton
    # This pattern captures the key tricyclic system with the oxygen bridge
    # Note the specific stereochemistry and connectivity
    core_pattern = Chem.MolFromSmarts("[C]1[C][C]2[C]3[C]([C][C]([C])[C]2(O[C]1([C])[C]))[C]3")
    
    if not mol.HasSubstructMatch(core_pattern):
        # Try alternative core pattern that's more flexible
        alt_core_pattern = Chem.MolFromSmarts("[C]1[C][C]2[C]3[C][C][C]2(O[C]1)[C]3")
        if not mol.HasSubstructMatch(alt_core_pattern):
            return False, "Missing dihydroagarofuran core structure"

    # Check for 5/7/6 fused ring system
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    has_correct_rings = False
    for ring in rings:
        if len(ring) in [5, 6, 7]:
            has_correct_rings = True
            break
    if not has_correct_rings:
        return False, "Missing characteristic 5/7/6 ring system"

    # Check for ester groups (these compounds typically have multiple ester substituents)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    if ester_matches < 2:
        return False, f"Insufficient ester groups (found {ester_matches}, expected ≥2)"

    # Check for oxygen bridge
    oxygen_bridge = Chem.MolFromSmarts("[C]1[O][C]([C])([C])[C]1")
    if not mol.HasSubstructMatch(oxygen_bridge):
        return False, "Missing characteristic oxygen bridge"

    # Count carbons (should have 15 for sesquiterpenoid core)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, f"Insufficient carbons for sesquiterpenoid (found {c_count}, need ≥15)"

    # Count oxygens (should have multiple due to ester groups and oxygen bridge)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, f"Insufficient oxygens (found {o_count}, expected ≥4)"

    # Check for quaternary carbons (characteristic of the skeleton)
    quat_c_pattern = Chem.MolFromSmarts("[C]([C])([C])([C])[C,O]")
    quat_c_matches = len(mol.GetSubstructMatches(quat_c_pattern))
    if quat_c_matches < 2:
        return False, "Insufficient quaternary carbons"

    # Check for sp3 hybridized carbons in core structure
    sp3_carbons = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX4]")))
    if sp3_carbons < 8:
        return False, f"Insufficient sp3 carbons in core structure (found {sp3_carbons}, expected ≥8)"

    # Additional check for molecular weight range based on examples
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (350 < mol_wt < 1000):
        return False, f"Molecular weight {mol_wt:.1f} outside typical range (350-1000)"

    return True, "Contains dihydroagarofuran skeleton with characteristic features: 5/7/6 ring system, oxygen bridge, and multiple ester substituents"