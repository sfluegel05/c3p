"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    A 3-oxo-fatty acyl-CoA is a fatty acyl-CoA structure with an oxo group on the third carbon.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a 3-oxo-fatty acyl-CoA, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for CoA structure, capturing nucleotide and phosphate linkage
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)C[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H](COP(=O)(O)[C@H]2O[C@H]([C@H](O)[C@H]2O)OP(=O)(O)O)n3cnc4c(N)ncnc43")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA structure found"

    # Look for oxo group on the third carbon of the fatty acyl chain
    oxo_fatty_acyl_pattern = Chem.MolFromSmarts("C[C@H](O)C(C)(C)COP(=O)(O)O")
    oxo_group = Chem.MolFromSmarts("O=CC(=O)C")
    if not mol.HasSubstructMatch(oxo_fatty_acyl_pattern) or not mol.HasSubstructMatch(oxo_group):
        return False, "No 3-oxo group found on the fatty acyl chain"

    return True, "Contains CoA structure and 3-oxo group on the fatty acyl chain"

# Example usage
example_smiles = "CC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)N1C=NC2=C1N=CN=C2N"
print(is_3_oxo_fatty_acyl_CoA(example_smiles))