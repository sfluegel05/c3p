"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    A 3-oxo-fatty acyl-CoA is characterized by a specific coenzyme A structure and
    an oxo group on the third carbon of the fatty acyl chain.

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

    # Define CoA structure using a specific SMARTS pattern
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)C[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC1COC(CO1)n1cnc2c(ncnc12)N")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA structure found"

    # Define 3-oxo group on the fatty acyl chain
    oxo_group_pattern = Chem.MolFromSmarts("CC(=O)C(=O)")
    if not mol.HasSubstructMatch(oxo_group_pattern):
        return False, "No 3-oxo group found on the fatty acyl chain"

    return True, "Contains CoA structure and 3-oxo group on the fatty acyl chain"