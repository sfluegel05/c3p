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
    
    # Update Coenzyme A pattern, ensuring it covers known cases accurately.
    coa_smarts = (
        "C(=O)SCCNC(=O)CCNC(=O)C[C@@H](O)C(C)COP(=O)(O)O"
        "P(=O)(O)OC[C@H]1O[C@H]([C@@H](O)[C@H](O)C1)n2cnc3c(N)ncnc23"
    )  # Expanding based on flexible structures in examples
    coa_pattern = Chem.MolFromSmarts(coa_smarts)

    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"
    
    # Define improved pattern for identifying branched carbon chains
    # Considering primary and secondary branching as a versatile descriptor.
    branched_smarts = "[CX3,CX4]([C])([C])[C]"  # Adjusting pattern to catch various branched structures
    branched_pattern = Chem.MolFromSmarts(branched_smarts)
    
    if not mol.HasSubstructMatch(branched_pattern):
        return False, "No branched-chain found in the fatty acid ligand"

    return True, "Contains both CoA moiety and appropriate branched-chain in fatty acid"