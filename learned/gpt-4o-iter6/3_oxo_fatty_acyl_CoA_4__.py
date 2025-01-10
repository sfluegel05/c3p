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

    # 3-oxo-fatty acyl pattern
    # Improved pattern to account for connectivity with CoA
    oxo_fatty_acyl_pattern = Chem.MolFromSmarts("[CX3](=O)[CX4][CX3](=O)SCC")
    if not mol.HasSubstructMatch(oxo_fatty_acyl_pattern):
        return False, "No connected 3-oxo-fatty acyl group found"

    # Enhanced CoA moiety pattern reflecting stereochemistry and more accurate structure
    # The pattern includes the CoA core with attached phosphates and stereocenters
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])([O-])=O.OP([O-])([O-])=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety with necessary stereochemistry not found"

    # Check for presence of 3 deprotonated phosphate groups
    deprotonated_phosphate_pattern = Chem.MolFromSmarts("P(O-)(O-)=O")
    phosphate_matches = mol.GetSubstructMatches(deprotonated_phosphate_pattern)
    if len(phosphate_matches) < 3:
        return False, "Insufficient deprotonated phosphate groups found"

    return True, "Contains 3-oxo-fatty acyl chain with CoA moiety including deprotonated phosphates"