"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    A 3-oxo-fatty acyl-CoA is a fatty acyl-CoA with an oxo group on the third carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the CoA pattern with detailed nucleotide and phosphate linkage
    coa_pattern = Chem.MolFromSmarts("C1(=O)COP(=O)(O)O[C@H]1O[C@H]2O[C@H](COP(=O)(O)O[C@H]3O[C@H]([C@H](O)[C@@H]3OP(=O)(O)O)n4cnc5c(N)ncnc54)C(O)C(O)C2OP(=O)(O)O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA structure found"

    # Look for 3-oxo linkage on the third carbon of a fatty acyl chain
    oxo_fatty_acyl_pattern = Chem.MolFromSmarts("C(=O)C(=O)SCCNC(=O)")
    if not mol.HasSubstructMatch(oxo_fatty_acyl_pattern):
        return False, "No 3-oxo-fatty acyl linkage found at the correct position"

    return True, "Contains CoA structure and 3-oxo group on the fatty acyl chain"

# Example usage
example_smiles = "CC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)N1C=NC2=C1N=CN=C2N"
print(is_3_oxo_fatty_acyl_CoA(example_smiles))