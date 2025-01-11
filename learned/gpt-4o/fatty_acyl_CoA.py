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
        bool: True if molecule is a fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string using RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for CoA moiety and thioester linkage
    # CoA pattern based on presence of adenine and phosphate groups
    coa_pattern = Chem.MolFromSmarts("n1cnc2c(n1)ncnc2N[CH2]CO[P]([O-])(=O)O[C@@H]1O[C@H]([C@@H](COP([O-])(=O)OP([O-])=O)[C@H]1O)O")
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)[O][C@H](O)C(C)(C)C")

    # Check for CoA structure
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing Coenzyme A moiety"

    # Check for thioester linkage
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Missing thioester linkage"
    
    # The presence of the CoA moiety and thioester bond confirms fatty acyl-CoA structure
    return True, "Contains Coenzyme A and thioester linkage to a fatty acid chain"