"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_long_chain_fatty_acyl_CoA_4_(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a long carbon chain typical of fatty acid backbone
    carbon_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCCC")  # At least 14 contiguous carbons
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No long carbon chain typical of fatty acid backbone found"
    
    # Look for CoA moiety pattern (thioester linkage -SCCNC(=O)...)
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety detected"
    
    # Check for presence of phosphate groups, indicative of the polyanionic state
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)([O-])")
    phosphate_count = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_count < 3:  # A typical CoA has multiple phosphate groups
        return False, f"Insufficient phosphate groups detected, found {phosphate_count}"
    
    return True, "Matches structure of long-chain fatty acyl-CoA(4-) with CoA moiety and long carbon chain"