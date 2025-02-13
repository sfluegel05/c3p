"""
Classifies: CHEBI:26244 prenols
"""
"""
Classifies: CHEBI:27290 prenol
A prenol is any alcohol possessing the general formula H-[CH2C(Me)=CHCH2]nOH, where the carbon skeleton
is composed of one or more isoprene units (biogenetic precursors of the isoprenoids).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prenol(smiles: str) -> tuple[bool, str]:
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

    # Look for alcohol group
    alcohol_pattern = Chem.MolFromSmarts("[OX1H]")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No alcohol group found"

    # Look for isoprene units (CH2=C(CH3)CH=CH2)
    isoprene_pattern = Chem.MolFromSmarts("[CH2]=[C@H]([CH3])[CH]=[CH2]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if not isoprene_matches:
        return False, "No isoprene units found"

    # Check for linear carbon skeleton
    skeleton_pattern = Chem.MolFromSmarts("[CH2]~[CH]~[CH2]~[CH]~[CH2]")
    skeleton_matches = mol.GetSubstructMatches(skeleton_pattern)
    if not skeleton_matches:
        return False, "Carbon skeleton is not linear"

    # Check for correct stereochemistry
    stereochem_pattern = Chem.MolFromSmarts("[CH2]=[C@H]([CH3])/[CH]=[CH/@H]")
    stereochem_matches = mol.GetSubstructMatches(stereochem_pattern)
    if not stereochem_matches:
        return False, "Incorrect stereochemistry"

    # Count isoprene units and carbon chain length
    n_isoprene = len(isoprene_matches)
    chain_length = len(skeleton_matches[0]) + 1
    if chain_length != n_isoprene * 5:
        return False, "Incorrect chain length for the number of isoprene units"

    return True, "Molecule is a prenol"