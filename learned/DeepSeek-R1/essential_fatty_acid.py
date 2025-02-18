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
    Essential fatty acids are polyunsaturated with a carboxylic acid group and all cis double bonds.

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

    # Count double bonds
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bonds < 2:
        return False, f"Found {double_bonds} double bonds, need at least 2"

    # Check all double bonds are cis (Z) or unspecified
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            stereo = bond.GetStereo()
            if stereo == Chem.BondStereo.STEREOE:
                return False, "Trans double bond present"

    return True, "Carboxylic acid with ≥2 cis double bonds"