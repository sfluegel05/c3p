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

    # Check for comprehensive Coenzyme A pattern (includes more of the structure)
    coa_comprehensive_pattern = Chem.MolFromSmarts("C(=O)NCCSC[CH2]NC(=O)[CH2]NC(=O)[C@@H](O)C(C)(C)OP([O-])(=O)OP([O-])(=O)OC[C@@H]1O[C@H]([C@H]1O)COS")
    if not mol.HasSubstructMatch(coa_comprehensive_pattern):
        return False, "Coenzyme A moiety not found"
    
    # Check for thioester linkage (C(=O)S)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"
    
    # Use a recursive SMARTS or dynamic range for medium acyl chain (6-12)
    medium_acyl_chain_pattern = Chem.MolFromSmarts("[CH2]~[CH2]~[CH2]~[CH2]~[CH2]", 6, 12)
    c_chain_matches = mol.GetSubstructMatches(medium_acyl_chain_pattern)
    if len(c_chain_matches) < 1:
        return False, "No medium-chain fatty acyl chain found"

    return True, "Contains medium-chain fatty acyl linked to Coenzyme A"