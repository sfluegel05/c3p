"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA(4-) based on its SMILES string.
    This class includes molecules with a 3-oxo group, a fatty acyl chain, and a CoA moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3-oxo-fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for the 3-oxo group pattern (C=O) at the 3-position of a carbon chain.
    oxo_pattern = Chem.MolFromSmarts("CC(=O)C")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "3-oxo group not found"
    
    # Check for a significant fatty acyl chain (a long continuous carbon chain)
    chain_pattern = Chem.MolFromSmarts("C" * 10)  # Expecting at least 10 continuous carbon atoms
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Fatty acyl chain too short"
    
    # Check for CoA signature including critical components (adenosine, phosphate, etc.)
    coa_pattern = Chem.MolFromSmarts("NC1=NC=CN=C1C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(C(O)C3OP(=O)(O)O)n4cnc5c(N)ncnc45)C(O)C2O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA structure not found"
    
    return True, "Contains 3-oxo group, adequate fatty acyl chain, and CoA moiety"