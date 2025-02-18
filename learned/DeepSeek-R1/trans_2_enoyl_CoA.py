"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
"""
Classifies: CHEBI:141322 trans-2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.
    A trans-2-enoyl-CoA has a trans double bond between carbons 2 and 3 of the acyl chain,
    connected via a thioester to coenzyme A.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is trans-2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for presence of thioester group (S-C(=O)-)
    thioester_pattern = Chem.MolFromSmarts("[SX2]-[CX3](=O)")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester group detected"

    # Check for CoA structural components (simplified to pantetheine and adenine parts)
    # Pattern matches S-C(=O)-C-N-C(=O) which is part of the CoA structure
    coa_pattern = Chem.MolFromSmarts("[SX2]-[CX3](=O)-[C]-[N]-[CX3](=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not detected"

    # Check for trans double bond at position 2,3 in acyl chain
    # SMARTS pattern matches C(=O)-C(/C=C/...)
    trans_pattern = Chem.MolFromSmarts("[CX3](=O)-[CX4H2]/[CX3H1]=[CX3H1]")
    if mol.HasSubstructMatch(trans_pattern):
        return True, "Trans-2-enoyl group connected to CoA via thioester"

    # Check alternative trans configuration pattern
    trans_pattern2 = Chem.MolFromSmarts("[CX3](=O)-[CX4H2]\\[CX3H1]=[CX3H1]")
    if mol.HasSubstructMatch(trans_pattern2):
        return True, "Trans-2-enoyl group connected to CoA via thioester"

    return False, "No trans-2-enoyl group found in acyl chain"