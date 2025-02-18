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

    # Improved CoA detection through multiple components
    # 1. Pantetheine chain pattern: S-C(=O)-C-C-N-C(=O)
    pantetheine_pattern = Chem.MolFromSmarts("[SX2]-[CX3](=O)-[CH2]-[CH2]-[NH]-[CX3](=O)")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Pantetheine chain not detected"
    
    # 2. Adenine detection from CoA structure
    adenine_pattern = Chem.MolFromSmarts("[n]1cnc2ncnc12")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Adenine moiety not found"

    # 3. Phosphate group verification
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX1-])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate groups not detected"

    # Check for trans double bond at position 2,3 in acyl chain
    # Find thioester carbonyl carbon
    thio_matches = mol.GetSubstructMatches(thioester_pattern)
    for match in thio_matches:
        s_idx, c_idx = match
        carbonyl_c = mol.GetAtomWithIdx(c_idx)
        # Get neighbors of carbonyl carbon (should be S and acyl chain carbon)
        for neighbor in carbonyl_c.GetNeighbors():
            if neighbor.GetIdx() != s_idx:  # Acyl chain carbon
                # Check next bond for trans double bond
                for bond in neighbor.GetBonds():
                    if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetBeginAtomIdx() == neighbor.GetIdx():
                        # Check stereochemistry
                        if bond.GetStereo() in [Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOTRANS]:
                            return True, "Trans double bond adjacent to thioester"
                        elif bond.GetStereo() in [Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOCIS]:
                            return False, "Cis double bond present instead of trans"

    # Fallback SMARTS pattern for trans configuration
    trans_pattern = Chem.MolFromSmarts("[CX3](=O)-[CH2]/[CX3]=[CX3]")
    if mol.HasSubstructMatch(trans_pattern):
        return True, "Trans-2-enoyl group detected via SMARTS"

    trans_pattern2 = Chem.MolFromSmarts("[CX3](=O)-[CH2]\\[CX3]=[CX3]")
    if mol.HasSubstructMatch(trans_pattern2):
        return True, "Trans-2-enoyl group detected via SMARTS"

    return False, "No trans-2-enoyl group found in acyl chain"