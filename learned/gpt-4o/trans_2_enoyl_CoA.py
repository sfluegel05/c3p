"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
from rdkit import Chem

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trans-2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for Coenzyme A backbone pattern
    coA_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coA_pattern):
        return False, "Coenzyme A moiety is missing or incorrect"
    
    # Check for trans-2-enoyl moiety pattern, specifically looking for the trans double bond 
    trans_2_enoyl_pattern = Chem.MolFromSmarts("C/C=C\\C(=O)")
    if trans_2_enoyl_pattern is None or not mol.HasSubstructMatch(trans_2_enoyl_pattern):
        return False, "No trans-2-enoyl moiety found"

    # Thioester linkage needs to be verified, i.e., an S-C(=O)-pattern
    thioester_pattern = Chem.MolFromSmarts("SC(=O)C")
    if thioester_pattern is None or not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    return True, "Detected Coenzyme A moiety with trans-2-enoyl moiety and appropriate thioester linkage"