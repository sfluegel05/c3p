"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA(4-) based on its SMILES string.
    This class includes molecules with a 3-oxo group, a long carbon chain typical for
    fatty acids, and a CoA moiety.

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
    
    # Refined 3-oxo pattern: might need flexibility around the adjacent carbon positions
    oxo_pattern = Chem.MolFromSmarts("C(=O)C")  # Just checking for the presence of the 3-oxo group
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "3-oxo group not found"
    
    # Improved fatty acyl chain pattern: more flexible pattern to accommodate longer chains
    chain_pattern = Chem.MolFromSmarts("CCCCCCCCC")  # 9 consecutive carbons indicating a long chain
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Fatty acyl chain too short"
    
    # Simplifying CoA structure match by identifying unique components of CoA
    # Example fragment: identifiable adenosine and phosphates
    coa_pattern = Chem.MolFromSmarts("NC=NC1C=CN=C1OC2COC(C2COP(=O)(O)OP(=O)(O)O)C(OP(=O)(O)O)n3cnc4c(N)ncnc34")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA structure not found"
    
    return True, "Contains 3-oxo group, adequate fatty acyl chain, and CoA moiety"