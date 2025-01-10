"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Check for 3-oxo-fatty acyl pattern - refined to include keto group
    oxo_fatty_acyl_pattern = Chem.MolFromSmarts("C(=O)C(=O)")
    if not mol.HasSubstructMatch(oxo_fatty_acyl_pattern):
        return False, "No 3-oxo-fatty acyl group found"

    # Improved check for CoA moiety specific pattern
    coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)OP([O-])=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Check for chirality and extended parts of CoA
    extended_pattern = Chem.MolFromSmarts("[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]n2cnc3c(N)ncnc32")
    if not mol.HasSubstructMatch(extended_pattern):
        return False, "Extended CoA parts with chirality not found"

    # Ensure presence of at least 2 phosphates with correct deprotonation
    phosphate_pattern = Chem.MolFromSmarts("OP([O-])(=O)[O-]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 2:
        return False, "Insufficient deprotonated phosphate groups found"

    return True, "Contains 3-oxo-fatty acyl chain with CoA moiety including deprotonated phosphates"