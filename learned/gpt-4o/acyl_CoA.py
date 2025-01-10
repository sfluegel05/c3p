"""
Classifies: CHEBI:17984 acyl-CoA
"""
from rdkit import Chem

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA is characterized by a -C(=O)S- linkage to a coenzyme A structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the thioester linkage (-C(=O)S-)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if thioester_pattern is None:
        return False, "Invalid thioester SMARTS pattern"
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found (-C(=O)S-)"

    # Define a SMARTS pattern for major components of the coenzyme A
    coa_pattern = Chem.MolFromSmarts("NCC(=O)CNC(=O)[C@H](O)C(C)(C)")
    if coa_pattern is None:
        return False, "Invalid CoA SMARTS pattern"
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No coenzyme A backbone found (missing typical CoA structure)"
    
    # Check for the terminal phosphate group combination commonly found in CoA
    phosphate_pattern = Chem.MolFromSmarts("COP(O)(=O)OP(O)(=O)O")
    if phosphate_pattern is None:
        return False, "Invalid phosphate SMARTS pattern"
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No terminal phosphate group pattern found"

    return True, "Molecule contains thioester linkage to a coenzyme A moiety"

# Example usage (not part of the final implementation):
# print(is_acyl_CoA("CCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"))