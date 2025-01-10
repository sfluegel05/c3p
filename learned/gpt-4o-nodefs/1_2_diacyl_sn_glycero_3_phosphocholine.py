"""
Classifies: CHEBI:57643 1,2-diacyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem

def is_1_2_diacyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 1,2-diacyl-sn-glycero-3-phosphocholine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is 1,2-diacyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol backbone pattern (C-C-C with oxygens)
    glycerol_pattern = Chem.MolFromSmarts("[CH2]C[CH2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for 2 ester groups attached to the glycerol backbone
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0][C]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Look for a phosphocholine group
    phosphocholine_pattern = Chem.MolFromSmarts("P(=O)(OCC[N+](C)(C)C)[O-]")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"

    # If all conditions are satisfied
    return True, "Contains glycerol backbone with 1,2-diacyl and phosphocholine group"

# Testing the function with provided SMILES strings
test_smiles = "C(C[N+](C)(C)C)OP(=O)([O-])OC[C@H](OC(CCCCCCCCCCCCCCC)=O)COC(=O)CCCCCCCCCCC"
print(is_1_2_diacyl_sn_glycero_3_phosphocholine(test_smiles))