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
        return False, "Invalid SMILES string"

    # 1. Check for phosphopantetheine moiety
    coA_pattern = Chem.MolFromSmarts("OC(C)C(NC(=O)CCNC(=O)S)COP(=O)([O-])OP(=O)([O-])")
    if not mol.HasSubstructMatch(coA_pattern):
        return False, "Missing phosphopantetheine and diphosphate connections"

    # 2. Possibly match adenine-related structure
    adenine_pattern = Chem.MolFromSmarts("N1C=NC2=C(N=CN=C2N1)")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Adenine-related structure not found"

    # 3. Check acyl chain with 3-position branching
    acyl_pattern = Chem.MolFromSmarts("CCC(=O)")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "Missing 3-substituted propionyl linkage"

    # 4. Check for negative charge states
    neg_charge_pat = Chem.MolFromSmarts("[O-]")
    neg_charge_matches = mol.GetSubstructMatches(neg_charge_pat)
    if len(neg_charge_matches) < 3:
        return False, "Insufficient deprotonated groups for 4- charge"

    return True, "Contains features of 3-substituted propionyl-CoA(4-)"