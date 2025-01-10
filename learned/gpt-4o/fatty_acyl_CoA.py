"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
from rdkit import Chem

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    A fatty acyl-CoA is defined by the presence of a CoA moiety and thioester bond linking to a fatty acid chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string using RDKit
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for CoA - broad search to accommodate structural variability
    coa_pattern = Chem.MolFromSmarts("n1c([nH]c2c(=O)ncnc2=O)c2c1ncn2C[C@H]1O[C@H](COP(=O)(O)O)[C@H](O)[C@H]1OP(=O)(O)O")
    
    # SMARTS pattern for thioester linkage
    thioester_pattern = Chem.MolFromSmarts("C(=O)SC")
    
    # Check for presence of the CoA structure
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing Coenzyme A moiety"

    # Check for presence of the thioester linkage
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Missing thioester linkage"
    
    return True, "Contains Coenzyme A moiety and thioester linkage to a fatty acid chain"