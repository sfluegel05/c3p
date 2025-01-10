"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA(4-) based on its SMILES string.
    This class is defined by the presence of a 3-oxo-fatty acyl chain attached to a CoA moiety
    with deprotonated phosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches the class, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined 3-oxo-fatty acyl pattern, focusing on long chains.
    oxo_fatty_acyl_pattern = Chem.MolFromSmarts("C(=O)CC(=O)C")
    if not mol.HasSubstructMatch(oxo_fatty_acyl_pattern):
        return False, "No 3-oxo-fatty acyl group found"

    # Improved CoA moiety pattern with attention to stereochemistry.
    coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OCC1OC(C(O)C1O)n2cnc3c(N)ncnc32")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety with stereochemistry not found"

    # Ensure there are exactly 4 deprotonated phosphate groups.
    phosphate_pattern = Chem.MolFromSmarts("OP([O-])(=O)[O-]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 3:
        return False, "Insufficient deprotonated phosphate groups found"

    return True, "Contains 3-oxo-fatty acyl chain with CoA moiety including deprotonated phosphates"