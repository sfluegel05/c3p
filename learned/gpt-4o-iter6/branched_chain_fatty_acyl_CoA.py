"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule belongs to the class of branched-chain fatty acyl-CoAs based on its SMILES string.
    A branched-chain fatty acyl-CoA has any branched-chain fatty acid attached to Coenzyme A via a thioester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule fits the class, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Refined Coenzyme A pattern covering more structural diversity
    coa_patterns = [
        "NC(=O)CCNC(=O)C[C@@H](O)C(C)(C)COP(=O)(O)O",  # Part of the molecule closer to base
        "n1cnc2c(N)ncnc12"  # adenine ring found at the very end
    ]
    # Check multiple parts of CoA to ensure full identification
    for pattern in coa_patterns:
        substruct = Chem.MolFromSmarts(pattern)
        if not mol.HasSubstructMatch(substruct):
            return False, f"No Coenzyme A part with pattern {pattern} found"

    # Improved pattern for identifying a branched carbon chain
    branched_smarts = "[CX3,CX4]([C])([C])[C]"  # general branched C
    branched_pattern = Chem.MolFromSmarts(branched_smarts)
    
    if not mol.HasSubstructMatch(branched_pattern):
        return False, "No branched-chain found in the fatty acid ligand"

    return True, "Contains both CoA moiety and appropriate branched-chain in fatty acid"