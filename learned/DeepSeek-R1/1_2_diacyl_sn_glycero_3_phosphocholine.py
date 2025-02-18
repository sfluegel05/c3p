"""
Classifies: CHEBI:57643 1,2-diacyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: 1,2-diacyl-sn-glycero-3-phosphocholine (conjugate base)
"""
from rdkit import Chem

def is_1_2_diacyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 1,2-diacyl-sn-glycero-3-phosphocholine based on its SMILES string.
    The structure must have:
    - sn-glycerol backbone with phosphocholine at position 3
    - Two esterified fatty acids at positions 1 and 2
    - Phosphate group in conjugate base form (negative charge)
    - Choline group (N+(C)(C)C)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # SMARTS pattern for core structure:
    # Chiral glycerol (sn-3) with two esters and phosphocholine
    # [C@] indicates stereochemistry matching sn-glycero-3 configuration
    pattern = Chem.MolFromSmarts(
        "[C@](COC(=O)*)(OC(=O)*)([H])OP(=[O])([O-])OCC[N+](C)(C)C"
    )
    if not mol.HasSubstructMatch(pattern):
        return False, "Core structure not matched"
    
    # Verify exactly two ester groups (already enforced by SMARTS)
    # Additional checks for fatty acid chains could be added here if needed
    
    return True, "1,2-diacyl-sn-glycero-3-phosphocholine structure"