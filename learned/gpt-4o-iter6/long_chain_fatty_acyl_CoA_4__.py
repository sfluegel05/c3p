"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_long_chain_fatty_acyl_CoA_4_(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.
    The compound should have a long-chain fatty acid linked to Coenzyme A via a thioester bond,
    with deprotonated phosphate groups.

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

    # Look for CoA thioester linkage pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Look for fatty acid chain (long carbon chain)
    long_chain_fatty_acid_pattern = Chem.MolFromSmarts("C" * 12 + "C(=O)SCC")
    if not mol.HasSubstructMatch(long_chain_fatty_acid_pattern):
        return False, "Insufficient carbon chain length"

    # Look for CoA backbone pattern including deprotonated phosphate groups
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA backbone or deprotonated phosphate groups not found"

    return True, "Matches long-chain fatty acyl-CoA(4-) structure with deprotonated phosphate groups"