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

    # Refined CoA pattern including pantetheine, diphosphate, and adenine
    coA_full_pattern = Chem.MolFromSmarts(
        "NC(=O)CCNC(=O)C[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H]([C@H](O)[C@H]1OP(=O)([O-])[O-])n2cnc3c(N)ncnc23"
    )
    if not mol.HasSubstructMatch(coA_full_pattern):
        return False, "Missing or incomplete coenzyme A structure"

    # 3-substituted propionyl chain with potential varying substitution
    propionyl_pattern = Chem.MolFromSmarts("SC(=O)CC")
    if not mol.HasSubstructMatch(propionyl_pattern):
        return False, "Missing or incorrectly structured 3-substituted propionyl chain"

    # Verify expected negative charges count
    # A propionyl-CoA(4-) should have exactly four negatively charged oxygens
    n_neg_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and atom.GetFormalCharge() == -1)
    if n_neg_oxygens != 4:
        return False, f"Incorrect number of negatively charged oxygens, found {n_neg_oxygens}, expected 4"

    return True, "Matches 3-substituted propionyl-CoA(4-) structure"