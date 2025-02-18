"""
Classifies: CHEBI:57643 1,2-diacyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: 1,2-diacyl-sn-glycero-3-phosphocholine (conjugate base)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_2_diacyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 1,2-diacyl-sn-glycero-3-phosphocholine.
    Requirements:
    - Correct sn-glycerol stereochemistry (C@@)
    - Two ester groups at positions 1 and 2
    - Phosphocholine group at position 3 with deprotonated phosphate
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Core structure with correct stereochemistry
    core_pattern = Chem.MolFromSmarts(
        "[C@@]([CH2]OC(=O)*)([CH2]OC(=O)*)([H])OP(=[O])([O-])OCC[N+](C)(C)C"
    )
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core structure not matched"
    
    # Verify two ester groups (O-C(=O)) attached to glycerol
    ester_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[O;X2]C(=O)"))
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"
    
    # Check fatty acid chains are present (longer than 2 carbons)
    # Example: each ester group should have at least 3 carbons (O-C-C-C=O)
    # This is a simplified check
    fatty_acid_pattern = Chem.MolFromSmarts("[O;X2][CX4][CX4][CX3]=O")
    if len(mol.GetSubstructMatches(fatty_acid_pattern)) < 2:
        return False, "Insufficient chain length for fatty acids"
    
    return True, "1,2-diacyl-sn-glycero-3-phosphocholine structure"