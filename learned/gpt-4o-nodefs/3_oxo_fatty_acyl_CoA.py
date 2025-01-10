"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    Key features:
    - Contains complete Coenzyme A moiety
    - Contains a thioester linkage
    - Contains a 3-oxo group in the fatty acyl chain

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # More detailed Coenzyme A moiety pattern including adenine, ribose, phosphates, and pantetheine
    coa_pattern = Chem.MolFromSmarts("CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H](COP(O)(=O)O)[C@@H](O)[C@@H]1OP(O)(O)=O)1[nH]c2c(N)ncnc2n1")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No full Coenzyme A moiety found"
    
    # Thioester linkage pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"
    
    # 3-oxo group pattern at third position of acyl chain
    oxo_pattern = Chem.MolFromSmarts("C(=O)CC(=O)")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No 3-oxo group at appropriate position found"
    
    return True, "Contains full CoA moiety, thioester linkage, and 3-oxo group in the fatty acyl chain"