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
    
    # Coenzyme A structural elements
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    ribose_pattern = Chem.MolFromSmarts("C1OC(C(O)COP(=O)(O)O)(OP(=O)(O)O)C(O)C1O")
    pantetheine_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@@H](O)C(C)(C)C")  # extended to include full structure
    
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine ring in CoA"
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "Missing ribose phosphate structure in CoA"
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing pantetheine structure in CoA"

    # Thioester linkage check
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCN")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage present"

    # Branched fatty acid chain
    branched_smarts = "[CX3;H2,CX3;H1]([CX4,CX3]([C])([C])[C])[CX4,CX3]([C])([C])[C]"  # more flexible branching pattern
    branched_pattern = Chem.MolFromSmarts(branched_smarts)
    
    if not mol.HasSubstructMatch(branched_pattern):
        return False, "No appropriate branched-chain found in the fatty acid"

    return True, "Contains both CoA moiety and branched-chain in fatty acid"