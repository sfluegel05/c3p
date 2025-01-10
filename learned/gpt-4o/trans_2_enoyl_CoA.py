"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Coenzyme A (CoA) structure pattern
    coA_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coA_pattern):
        return False, "No Coenzyme A moiety found"

    # Trans-2-enoyl moiety pattern
    trans_2_enoyl_pattern = Chem.MolFromSmarts("C\C=C\(C(=O)")
    if not mol.HasSubstructMatch(trans_2_enoyl_pattern):
        return False, "No trans-2-enoyl moiety found"

    # Thioester linkage (S-C(=O)-C)
    thioester_pattern = Chem.MolFromSmarts("S[C](=O)")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) == 0:
        return False, "No thioester linkage found"

    return True, "Contains Coenzyme A moiety with trans-2-enoyl moiety and thioester linkage, qualifying as trans-2-enoyl-CoA"