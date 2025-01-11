"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for thioester linkage (C(=O)S)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"
    
    # Check for Coenzyme A moiety
    coa_pattern = Chem.MolFromSmarts("NCC(=O)CC(=O)[C@H](O)C(C)(C)COP(=O)(O)OCC1OC(C(O)C1OP(=O)(O)O)n2cnc3c(N)cnc(N)c23")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"
    
    # Check for medium-chain fatty acid length (6 to 12 carbons)
    carbon_chain_pattern = Chem.MolFromSmarts("[#6]~[#6]~[#6]~[#6]~[#6]~[#6,7,#6,#6,#6,#6,#6,#6,#6,#6]") # 6 to 12 carbons
    carbon_chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    if not carbon_chain_matches:
        return False, "Aliphatic chain length not within medium-chain range (6-12 carbons)"

    return True, "Molecule is a medium-chain fatty acyl-CoA with CoA moiety and appropriate chain length"

# Example SMILES for testing:
example_smiles = "CCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
print(is_medium_chain_fatty_acyl_CoA(example_smiles))