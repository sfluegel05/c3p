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
from rdkit.Chem import rdMolDescriptors

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
    chroman_pattern = Chem.MolFromSmarts("[OX2;R]1[CX4;R]2=[CX3]([CH3])[CX3]=[CX3][CX3]=[CX3]2[CX4]=[CX3][CX3]1")
    if not mol.HasSubstructMatch(chroman_pattern):
        return False, "No chroman-6-ol skeleton found"

    # Look for hydrocarbon isoprenoid chain at position 2
    isoprenoid_pattern = Chem.MolFromSmarts("[CX4]([CX4]([CX4]=[CX3][CX3]=[CX3])[CX4]=[CX3][CX3]=[CX3])[CX4]=[CX3][CX3]=[CX3]")
    isoprenoid_match = mol.GetSubstructMatch(isoprenoid_pattern)
    if not isoprenoid_match:
        return False, "No isoprenoid chain at position 2"

    # Check if chain is saturated or triply unsaturated
    unsaturated_bonds = sum(1 for bond in mol.GetBondBetweenAtoms(isoprenoid_match[0], isoprenoid_match[1]).GetBondTypeAsDouble())
    if unsaturated_bonds not in [0, 3]:
        return False, "Side chain should be either saturated or triply unsaturated"

    # Check length of isoprenoid chain
    chain_length = 0
    for i in range(len(isoprenoid_match) - 1):
        atom1 = mol.GetAtomWithIdx(isoprenoid_match[i])
        atom2 = mol.GetAtomWithIdx(isoprenoid_match[i + 1])
        if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
            chain_length += 1
    if chain_length != 3:
        return False, f"Isoprenoid chain length is {chain_length}, should be 3"

    # Check number of carbon and hydrogen atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    if c_count < 20 or c_count > 35 or h_count < 30 or h_count > 55:
        return False, "Number of carbon and hydrogen atoms outside expected range for tocols"

    return True, "Contains a chroman-6-ol skeleton substituted at position 2 with a saturated or triply-unsaturated 3-isoprenoid hydrocarbon chain"