"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: essential fatty acid (CHEBI:67011)
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    Essential fatty acids are polyunsaturated with a carboxylic acid group, all cis double bonds,
    and no other functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an essential fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"

    # Check for carboxylic acid group
    carboxyl = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl):
        return False, "No carboxylic acid group"

    # Check oxygen count is exactly 2 (from COOH)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 2:
        return False, f"Expected 2 oxygen atoms, found {o_count}"

    # Check all other atoms are carbon or hydrogen
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [6, 1, 8]:  # Allow oxygens in COOH
            return False, "Contains non-carbon/hydrogen atoms besides COOH"

    # Count all double bonds (including COOH)
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bonds < 2:
        return False, f"Found {double_bonds} double bonds, need at least 2"

    # Check all double bonds (except COOH) are cis or unspecified
    carboxyl_c = [match[0] for match in mol.GetSubstructMatches(carboxyl)]
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            # Skip COOH double bond
            atoms = [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]
            if any(idx in carboxyl_c for idx in atoms):
                continue
            # Check cis configuration
            if bond.GetStereo() == Chem.BondStereo.STEREOE:
                return False, "Trans double bond present"

    return True, "Polyunsaturated fatty acid with all cis double bonds and no other functional groups"