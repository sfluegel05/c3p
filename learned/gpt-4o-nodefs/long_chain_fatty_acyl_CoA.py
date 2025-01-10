"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    
    # Define the Coenzymes A structure pattern
    coa_pattern_smarts = "NC(=O)CCNC(=O)[C@H](O)C(C)COP(=O)(O)OP(=O)(O)O"
    coa_pattern = Chem.MolFromSmarts(coa_pattern_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing Coenzyme A structure"
        
    # Define the thiol ester connectivity pattern
    thiol_ester_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)")
    if not mol.HasSubstructMatch(thiol_ester_pattern):
        return False, "Missing proper thiol ester linkage"
    
    # Define a long hydrocarbon chain pattern, allowing for flexibility in length (16+ carbons)
    long_chain_pattern_smarts = "C(=O)C{14,}CCC"
    long_chain_pattern = Chem.MolFromSmarts(long_chain_pattern_smarts)
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No sufficiently long hydrocarbon chain found"

    # Check for unsaturations (double bonds)
    # This should be more flexible to allow the range of possibilities seen in sample SMILES.
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_count = len(mol.GetSubstructMatches(double_bond_pattern))
    if double_bond_count > 6:
        return False, "Too many double bonds for typical long-chain fatty acyl-CoA"
    
    return True, "Matches long-chain fatty acyl-CoA features"