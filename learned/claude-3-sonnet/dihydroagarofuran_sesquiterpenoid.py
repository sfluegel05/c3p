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
    # and specific connectivity characteristic of dihydroagarofuran
    core_pattern = Chem.MolFromSmarts("[C]1[C][C]2[C]3[C]([C][C]([C])[C]3(O[C]1([C])[C])[C]2)([C,O])")
    
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing dihydroagarofuran core structure"

    # Count rings to verify tricyclic system
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Insufficient number of rings"

    # Check for characteristic oxygen bridge
    oxygen_bridge = Chem.MolFromSmarts("[C]1[O][C]2[C][C]1[C]2")
    if not mol.HasSubstructMatch(oxygen_bridge):
        return False, "Missing characteristic oxygen bridge"

    # Check for ester groups (these compounds typically have multiple ester substituents)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    if ester_matches < 3:
        return False, f"Insufficient ester groups (found {ester_matches}, expected ≥3)"

    # Count carbons (should have at least 15 for sesquiterpenoid core)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, f"Insufficient carbons for sesquiterpenoid (found {c_count}, need ≥15)"

    # Count oxygens (should have multiple due to ester groups and oxygen bridge)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, f"Insufficient oxygens (found {o_count}, expected ≥6)"

    # Check for quaternary carbons (characteristic of the skeleton)
    quat_c_pattern = Chem.MolFromSmarts("[C]([C])([C])([C])[C]")
    quat_c_matches = len(mol.GetSubstructMatches(quat_c_pattern))
    if quat_c_matches < 2:
        return False, "Insufficient quaternary carbons"

    # Molecular weight check (based on example structures)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (400 < mol_wt < 900):
        return False, f"Molecular weight {mol_wt:.1f} outside typical range (400-900)"

    # Additional check for sp3 carbons in the core
    sp3_c = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX4]")))
    if sp3_c < 10:
        return False, f"Insufficient sp3 carbons (found {sp3_c}, expected ≥10)"

    return True, "Contains dihydroagarofuran skeleton with characteristic features: tricyclic system, oxygen bridge, and multiple ester substituents"