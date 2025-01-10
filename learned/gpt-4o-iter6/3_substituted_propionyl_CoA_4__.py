"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
from rdkit import Chem

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-substituted propionyl-CoA(4-) 
    based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-substituted propionyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string, cannot parse"

    # 1. Check for coenzyme A (CoA) backbone structure
    # This includes pantetheine connected via its phosphate groups
    coA_pattern = Chem.MolFromSmarts("C(C)(C)COP(=O)([O-])OP(=O)([O-])[O]c1ncnc2ncnc12")
    if not mol.HasSubstructMatch(coA_pattern):
        return False, "Missing coenzyme A structure"

    # 2. Check for 3-substituted propionyl chain
    # The 3-position should have the branching or substitution
    # This is a simplified assumption
    propionyl_pattern = Chem.MolFromSmarts("SC(=O)CCC")
    if not mol.HasSubstructMatch(propionyl_pattern):
        return False, "Missing 3-substituted propionyl chain"

    # 3. Check for correct negative charges
    # Verify the number of deprotonated phosphate oxygens
    neg_charge_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[O-]"))
    if len(neg_charge_matches) < 4:
        return False, f"Insufficient negative charges, found {len(neg_charge_matches)}"

    return True, "Matches 3-substituted propionyl-CoA(4-) structure"