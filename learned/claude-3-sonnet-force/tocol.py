"""
Classifies: CHEBI:39437 tocol
"""
"""
Classifies: CHEBI:53685 tocol

A tocol is a chromanol with a chroman-6-ol skeleton that is substituted at position 2 by a
saturated or triply-unsaturated hydrocarbon chain consisting of three isoprenoid units.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tocol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for chroman-6-ol skeleton
    chroman_pattern = Chem.MolFromSmarts("[OX2;R]1[CX4;R]2=[CX3]([CX3]=[CX3][CX3]=[CX3]2)[CX4]=[CX3][CX3]1")
    if not mol.HasSubstructMatch(chroman_pattern):
        return False, "No chroman-6-ol skeleton found"

    # Look for isoprenoid chain at position 2
    isoprenoid_pattern = Chem.MolFromSmarts("[CX4]([CX4]([CX4]([CX4]=[CX3][CX3]=[CX3])=[CX3][CX3]=[CX3])=[CX3][CX3]=[CX3])[CX4]=[CX3][CX3]=[CX3]")
    isoprenoid_match = mol.GetSubstructMatches(isoprenoid_pattern)
    if not isoprenoid_match:
        return False, "No isoprenoid chain at position 2"

    # Check if chain is saturated or triply unsaturated
    for match in isoprenoid_match:
        unsaturated_bonds = sum(1 for bond in mol.GetBondBetweenAtoms(match[0], match[1]).GetBondTypeAsDouble())
        if unsaturated_bonds not in [0, 3]:
            return False, "Side chain should be either saturated or triply unsaturated"

    return True, "Contains a chroman-6-ol skeleton substituted at position 2 with a saturated or triply-unsaturated 3-isoprenoid hydrocarbon chain"