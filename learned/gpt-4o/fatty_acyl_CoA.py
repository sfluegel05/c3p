"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
from rdkit import Chem

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    Fatty acyl-CoA is an acyl-CoA derived from the condensation of
    coenzyme A with a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for the thioester linkage (C(=O)S) between acyl group and CoA
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Adjusted pattern for a fatty acyl chain, accounting for flexible chain length and unsaturations
    fatty_acyl_pattern = Chem.MolFromSmarts("C(=O)[C,c;!R][C,c;!R][C,c;!R]")  # Allow more flexibility and unsaturation
    if not mol.HasSubstructMatch(fatty_acyl_pattern):
        return False, "No suitable fatty acyl chain found"

    # Improved coenzyme A structure SMARTS pattern, focusing on specific features
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)([O-])=O)n2cnc3nc(N)nc2-3")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA structure found"

    return True, "Contains thioester linkage with fatty acyl chain and CoA structure"