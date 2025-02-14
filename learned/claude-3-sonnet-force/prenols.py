"""
Classifies: CHEBI:26244 prenols
"""
"""
Classifies: CHEBI:36360 prenols

Prenols are defined as any alcohol possessing the general formula H-[CH2C(Me)=CHCH2]nOH, where the
carbon skeleton is composed of one or more isoprene units (biogenetic precursors of the isoprenoids).
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors


def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for isoprene units (C=CC(C)C) in both cis and trans configurations
    isoprene_pattern = Chem.MolFromSmarts("[CH2]=[C@H][CH2][CX4]([CH3])[CH2]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    isoprene_pattern_trans = Chem.MolFromSmarts("[CH2]=[CH][CH2][CX4]([CH3])[CH2]")
    isoprene_matches_trans = mol.GetSubstructMatches(isoprene_pattern_trans)
    if not isoprene_matches and not isoprene_matches_trans:
        return False, "No isoprene units found"

    # Look for hydroxyl group
    hydroxy_pattern = Chem.MolFromSmarts("[OX1H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if not hydroxy_matches:
        return False, "No hydroxyl group found"

    # Check for linear or branched carbon skeleton
    skeleton_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]")
    skeleton_matches = mol.GetSubstructMatches(skeleton_pattern)
    if not skeleton_matches:
        return False, "Carbon skeleton is not linear or branched"

    # Check if isoprene units are part of the carbon skeleton
    for match in isoprene_matches + isoprene_matches_trans:
        isoprene_atom_ids = [mol.GetAtomWithIdx(idx).GetIdx() for idx in match]
        if not any(set(isoprene_atom_ids).issubset(set(match)) for match in skeleton_matches):
            return False, "Isoprene units not part of the carbon skeleton"

    return True, "Contains isoprene units as part of a linear or branched carbon skeleton with a terminal hydroxyl group"