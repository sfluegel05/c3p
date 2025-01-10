"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
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

    # Define a more precise thioester linkage pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No CoA thioester linkage found"

    # Check for a long-chain fatty acid, minimum 12 contiguous carbons
    long_chain_fatty_acid_pattern = Chem.MolFromSmarts("[C]" + "{8,}" + "C(=O)S")
    if not mol.HasSubstructMatch(long_chain_fatty_acid_pattern):
        return False, "Insufficient carbon chain length in the fatty acid"

    # Look for CoA backbone with deprotonated phosphate groups
    coa_diphosphate_pattern = Chem.MolFromSmarts("COP([O-])(=O)OP([O-])(=O)O[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_diphosphate_pattern):
        return False, "CoA backbone or deprotonated phosphate groups not found"

    return True, "Matches long-chain fatty acyl-CoA(4-) structure with deprotonated phosphate groups"

# Example usage:
smiles = "CC(C)(COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)CCCCCCC\C=C/C\C=C/CC=C"
result, reason = is_long_chain_fatty_acyl_CoA_4__(smiles)
print(result, reason)