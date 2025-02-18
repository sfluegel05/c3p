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

    # Improved CoA pattern: S-C(=O)-C-C-N-C(=O) (part of pantetheine)
    coa_pattern = Chem.MolFromSmarts("[SX2]-[CX3](=O)-[CX4H2]-[CX4H2]-[NX3H1]-[CX3](=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not detected"

    # Check for trans double bond at position 2,3 in acyl chain (adjacent to thioester)
    # SMARTS pattern: thioester carbonyl followed by trans double bond
    # Pattern matches [S]-C(=O)-C=C with trans configuration
    trans_pattern = Chem.MolFromSmarts("[CX3](=O)-[CX4H2]/[CX3H1]=[CX3H1]")
    if mol.HasSubstructMatch(trans_pattern):
        return True, "Trans-2-enoyl group connected to CoA via thioester"

    trans_pattern2 = Chem.MolFromSmarts("[CX3](=O)-[CX4H2]\\[CX3H1]=[CX3H1]")
    if mol.HasSubstructMatch(trans_pattern2):
        return True, "Trans-2-enoyl group connected to CoA via thioester"

    # Alternative approach: Check all double bonds adjacent to thioester's carbonyl for trans configuration
    # Find the thioester carbonyl carbon
    thio_matches = mol.GetSubstructMatches(thioester_pattern)
    for match in thio_matches:
        s_idx, c_idx = match
        carbonyl_c = mol.GetAtomWithIdx(c_idx)
        # Get neighbors of carbonyl carbon (should be S and next carbon)
        for neighbor in carbonyl_c.GetNeighbors():
            if neighbor.GetIdx() != s_idx:  # This is the acyl chain carbon
                acyl_start = neighbor
                # Check if this carbon has a double bond
                for bond in acyl_start.GetBonds():
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        # Check if double bond is trans
                        if bond.GetStereo() == Chem.BondStereo.STEREOE:
                            return True, "Trans double bond adjacent to thioester"
                        elif bond.GetStereo() == Chem.BondStereo.STEREOZ:
                            return False, "Cis double bond adjacent to thioester"

    return False, "No trans-2-enoyl group found in acyl chain"