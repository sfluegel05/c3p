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

    # Define SMARTS for the thioester linkage (-C(=O)S-)
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCN")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found (-C(=O)S-)"

    # Define a comprehensive SMARTS pattern for coenzyme A structure
    coa_pattern = Chem.MolFromSmarts("[nH]1cnc2c1ncnc2N[C@H]3C(C)(C)O[C@@H](COP(O)(=O)O)O[C@H]3P(=O)(O)OP(=O)(O)O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No coenzyme A backbone found (missing typical CoA structure)"

    return True, "Molecule contains thioester linkage to a coenzyme A moiety"