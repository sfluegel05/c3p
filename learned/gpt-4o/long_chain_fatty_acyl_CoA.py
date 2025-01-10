"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for CoA moiety
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A (CoA) moiety not found"
    
    # Define SMARTS pattern for the thioester linkage indicating fatty acyl-CoA
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage required for acyl-CoA not found"
    
    # Find the fatty acyl component and check carbon chain length (13 to 22)
    # Use SMARTS to find a long chain of carbon atoms
    carbon_chain_pattern = Chem.MolFromSmarts("C{13,22}")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Fatty acyl chain not in the correct length range (13-22 carbons)"
    
    return True, "Structure matches long-chain fatty acyl-CoA"