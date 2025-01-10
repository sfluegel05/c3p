"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Classifies if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a broader CoA pattern
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@@H](O)C(C)(C)COP(=O)(O)OC[C@H]1O[C@H](O)[C@@H](OP(=O)(O)O1)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A structure not found"
    
    # Pattern for the thioester bond linking CoA
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester bond not found"
    
    # Pattern for branching: allow more generic branched structures
    branch_pattern = Chem.MolFromSmarts("C(C)(C)[C;!H1,XR]")
    if not mol.HasSubstructMatch(branch_pattern):
        return False, "No branched alkyl chain found"

    return True, "Contains branched-chain fatty acyl-CoA structure"