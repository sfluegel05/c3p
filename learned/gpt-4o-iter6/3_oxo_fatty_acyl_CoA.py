"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for 3-oxo-fatty acid (a ketone at position 3 on a carbon chain)
    oxo_fatty_acid_pattern = Chem.MolFromSmarts("C(=O)CC")
    if not mol.HasSubstructMatch(oxo_fatty_acid_pattern):
        return False, "No 3-oxo fatty acid group found"
    
    # Pattern for CoA moiety (capturing phosphopantetheine and 3'-phosphoadenosine diphosphate)
    coa_pattern = Chem.MolFromSmarts("[C@H](O)C(C)(C)COP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Look for the thioester linkage (C(=O)S-C, ester linkage with sulfur)
    thioester_pattern = Chem.MolFromSmarts("C(=O)SC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found between 3-oxo and CoA"

    return True, "Structure matches 3-oxo-fatty acyl-CoA with 3-oxo group, CoA moiety, and thioester linkage"