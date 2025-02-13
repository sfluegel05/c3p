"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.
    A medium-chain fatty acyl-CoA is a fatty acyl-CoA with a fatty acid chain of 6-12 C atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for CoA adenine part: a common fragment of Coenzyme A
    coa_adenine_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc2n1)C3CC(O[C@@H]4[C@H](O)C[C@@H](COP(O)(O)=O)O4)C3")
    if not mol.HasSubstructMatch(coa_adenine_pattern):
        return False, "Missing Coenzyme A adenine component"
    
    # Check for CoA phosphate component 
    coa_phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O[C@H]1[C@H](O)[C@@H](O)C(O)C1")
    if not mol.HasSubstructMatch(coa_phosphate_pattern):
        return False, "Missing Coenzyme A phosphate component"

    # Check for thioester linkage (should be part of the fatty acid-CoA linkage): -C(=O)-SCC
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)C")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found to CoA"

    # Check for medium-chain fatty acid length (6-12 carbons)
    # Adjust the carbon chain pattern to belong solely to the fatty acid component
    # The pattern denotes the carbon chain length following the C(=O) moiety
    chain_pattern_6_12 = Chem.MolFromSmarts("CCCCCC")
    carbon_chains = mol.GetSubstructMatches(chain_pattern_6_12)
    if any(6 <= len(chain) <= 12 for chain in carbon_chains):
        return True, "Contains medium-chain fatty acyl-CoA structure with appropriate length"
    
    return False, "Fatty acid chain not in medium range (6-12 carbons)"