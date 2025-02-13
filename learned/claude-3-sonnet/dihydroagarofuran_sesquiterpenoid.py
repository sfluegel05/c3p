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
    Dihydroagarofuran sesquiterpenoids have a characteristic tricyclic (5/6/6) ring system
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
    # This captures the key tricyclic system with specific connectivity
    core_pattern = Chem.MolFromSmarts("""
        [C]1[CH2][C]2[C]3[C]([C][C]([CH3])[C]2(O[C]1([CH3])[CH3]))[C]3
    """.replace(" ", ""))
    
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing dihydroagarofuran core structure"

    # Check for characteristic oxygen bridge in correct position
    oxygen_bridge = Chem.MolFromSmarts("[C]1[O][C]([CH3])([CH3])[C]1")
    if not mol.HasSubstructMatch(oxygen_bridge):
        return False, "Missing characteristic oxygen bridge"

    # Check for multiple ester groups (these compounds typically have 3 or more ester substituents)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    if ester_matches < 3:
        return False, f"Insufficient ester groups (found {ester_matches}, expected ≥3)"

    # Count carbons (should have 15 for sesquiterpenoid core plus additional for substituents)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, f"Insufficient carbons for sesquiterpenoid (found {c_count}, need ≥15)"

    # Count oxygens (should have multiple due to ester groups and oxygen bridge)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:  # Increased from 4 based on examples
        return False, f"Insufficient oxygens (found {o_count}, expected ≥6)"

    # Check for quaternary carbons (characteristic of the skeleton)
    quat_c_pattern = Chem.MolFromSmarts("[C]([C])([C])([C])[C,O]")
    quat_c_matches = len(mol.GetSubstructMatches(quat_c_pattern))
    if quat_c_matches < 3:  # Increased from 2 based on examples
        return False, "Insufficient quaternary carbons"

    # Check for specific ring sizes (5/6/6 system)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    ring_sizes = [len(ring) for ring in rings]
    if not (5 in ring_sizes and ring_sizes.count(6) >= 2):
        return False, "Missing characteristic 5/6/6 ring system"

    # Additional check for sp3 carbons in core structure
    sp3_pattern = Chem.MolFromSmarts("[CX4]")
    sp3_matches = len(mol.GetSubstructMatches(sp3_pattern))
    if sp3_matches < 10:  # Increased based on core structure requirements
        return False, f"Insufficient sp3 carbons (found {sp3_matches}, expected ≥10)"

    # Check molecular weight range based on examples
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (400 < mol_wt < 1000):  # Adjusted range based on examples
        return False, f"Molecular weight {mol_wt:.1f} outside typical range (400-1000)"

    return True, "Contains dihydroagarofuran skeleton with characteristic features: 5/6/6 ring system, oxygen bridge, and multiple ester substituents"