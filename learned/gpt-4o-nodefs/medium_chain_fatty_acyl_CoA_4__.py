"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improved Coenzyme A SMARTS pattern (more comprehensive)
    # Includes ribose, phosphate groups, pantetheine, and adenine moiety
    coa_comprehensive_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)CCSC[CH2]NC(=O)[CH2]OP([O-])(=O)OC[C@H]1O[C@H](n2cnc3c(N)ncnc32)[C@H]1O")
    if not mol.HasSubstructMatch(coa_comprehensive_pattern):
        return False, "Coenzyme A moiety not found"

    # Check for thioester linkage (C(=O)S)
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"
    
    # Check for medium acyl chain length (6 to 12 carbons)
    # This can be complex since acyl groups can be present in different conformations and forms
    # Examine the count in carbon chains linked via thioester pattern
    carbon_chain_pattern = Chem.MolFromSmarts("C([C;R0])*")
    chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    for match in chain_matches:
        chain_length = sum(1 for idx in match if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if 6 <= chain_length <= 12:
            return True, "Contains medium-chain fatty acyl linked to Coenzyme A"
    
    return False, "No medium-chain fatty acyl chain found"