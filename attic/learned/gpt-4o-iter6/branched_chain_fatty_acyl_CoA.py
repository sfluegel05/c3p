"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    A branched-chain fatty acyl-CoA has a branched-chain fatty acid attached to Coenzyme A via a thioester bond.

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

    # Refined Coenzyme A SMARTS pattern
    # Includes the thioester linkage, pantetheine phosphate, and adenine nucleotide portions as distinguished in CoA
    coa_smarts = (
        "C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O"
        "[C@H]([C@H](O)[C@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    )
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"

    # Pattern for branched fatty acids (at least one tertiary carbon)
    # A carbon with three non-hydrogen connected atoms (e.g., -C(CH3)3)
    branched_pattern = Chem.MolFromSmarts("[CX4](C)(C)C")
    
    if not mol.HasSubstructMatch(branched_pattern):
        return False, "No branched-chain in the fatty acid"

    return True, "Contains CoA moiety and branched-chain in the fatty acid"