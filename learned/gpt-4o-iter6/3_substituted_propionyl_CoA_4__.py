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

    # 1. Check for coenzyme A (CoA) elements, more generally defined.
    # This includes pantetheine connected via phosphate and adenine elements
    coA_pattern = Chem.MolFromSmarts("C(C)(C)COP(=O)([O-])O[P](=O)([O-])O[C@@H]1O[C@H]([C@@H](O)[C@H]1OP(=O)([O-])[O-])n2cnc3c(N)ncnc23")
    if not mol.HasSubstructMatch(coA_pattern):
        return False, "Missing or incomplete coenzyme A structure"

    # 2. Check for a 3-substituted propionyl chain
    # Match SC(=O)CCC* using a wildcard for the substitution at the third position
    propionyl_pattern = Chem.MolFromSmarts("SC(=O)CC[CX4]")
    if not mol.HasSubstructMatch(propionyl_pattern):
        return False, "Missing or incorrectly structured 3-substituted propionyl chain"

    # 3. Revalidate correct negative charges
    # Verify the number of deprotonated phosphate oxygens
    neg_charge_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[O-]"))
    if len(neg_charge_matches) < 4:
        return False, f"Insufficient negative charges, found {len(neg_charge_matches)}"

    return True, "Matches 3-substituted propionyl-CoA(4-) structure"