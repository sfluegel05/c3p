"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    This class is characterized by the formal condensation of the thiol group of coenzyme A
    with the carboxy group of any branched-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify the thioester linkage (C(=O)SCC) - this is part of typical acyl-CoA structures
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found indicative of acyl-CoA"

    # Identify key components of the Coenzyme A
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)N")
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)OP(=O)(O)O")
    if not mol.HasSubstructMatch(adenine_pattern) or not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Missing adenine or phosphate groups indicative of coenzyme A"

    # Look for branching in the fatty acyl chain
    # Check for tertiary or quaternary carbons
    branching_patterns = [
        Chem.MolFromSmarts("[C](C)(C)C"),
        Chem.MolFromSmarts("[C](C)(C)(C)")  # This accounts for more complex branching
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in branching_patterns):
        return False, "No branching in the fatty acid chain found"

    return True, "Structure matches branched-chain fatty acyl-CoA"