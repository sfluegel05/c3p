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
    
    # Define and detect Coenzyme A pattern using SMARTS
    coa_smarts = (
        "C(=O)SCCNC(=O)CCNC(=O)C[C@H](O)C(C)COP(O)(=O)OP(O)(=O)OCC1OC"
        "(n2cnc3c(ncnc23)N)[C@H]1O"
    )
    coa_pattern = Chem.MolFromSmarts(coa_smarts)

    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"
    
    # Define pattern for branched chains in fatty acids
    # A carbon atom connected to three other carbons
    branched_smarts = "[CX4](C)(C)C"  # Detects _branched_ carbon patterns: tertiary or quaternary carbons
    branched_pattern = Chem.MolFromSmarts(branched_smarts)
    
    if not mol.HasSubstructMatch(branched_pattern):
        return False, "No branched-chain in the fatty acid"

    return True, "Contains CoA moiety and branched-chain in the fatty acid"