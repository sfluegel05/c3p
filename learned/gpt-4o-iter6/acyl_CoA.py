"""
Classifies: CHEBI:17984 acyl-CoA
"""
from rdkit import Chem

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA is a thioester formed from the condensation of a carboxylic acid and coenzyme A.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for thioester bond pattern S-C(=O) as basic indicator
    thioester_pattern = Chem.MolFromSmarts("S(=O)[CX3]")  # More correctly specifies sulfur bonded to a carbonyl carbon
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"
    
    # Simplified pattern for major structural elements of Coenzyme A (e.g., pantetheine connected to adenine moiety)
    coa_pattern = Chem.MolFromSmarts("NCC(=O)CCNC(=O)[C@H](O)[C@H](C)[C@H](O)COP(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"

    return True, "Contains thioester linkage and Coenzyme A structure"