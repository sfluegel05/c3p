"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
from rdkit import Chem

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Classifies a molecule as trans-2-enoyl-CoA based on its SMILES string.
    Checks for a trans double bond, thioester linkage, and Coenzyme A moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as trans-2-enoyl-CoA
        str: Reason for the classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for confirming a trans double bond with plausible adjacent groups
    trans_double_bond_pattern = Chem.MolFromSmarts("C/C=C\\C(=O)S")
    if not mol.HasSubstructMatch(trans_double_bond_pattern):
        return False, "No suitable trans double bond with thioester found"
    
    # SMARTS for Coenzyme A motif
    coa_pattern = Chem.MolFromSmarts(
        "SC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)O" +
        "C[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(=O)O)n2cnc3c(N)ncnc23"
    )
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A component not detected correctly"
    
    return True, "Matches key features of trans-2-enoyl-CoA"