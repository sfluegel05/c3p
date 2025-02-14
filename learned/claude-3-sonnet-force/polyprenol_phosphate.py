"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
"""
Classifies: CHEBI:37412 polyprenol phosphate

A prenol phosphate resulting from the formal condensation of the terminal allylic 
hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for phosphate group (P(O)(O)=O)
    phosphate_pattern = Chem.MolFromSmarts("P(O)(O)=O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Look for linear isoprenoid chain (C=C-C=C-C=C...)
    isoprenoid_pattern = Chem.MolFromSmarts("C=CC(C)C=CC=CC=C")
    isoprenoid_matches = mol.GetSubstructMatches(isoprenoid_pattern)

    if not isoprenoid_matches:
        return False, "No linear isoprenoid chain found"

    # Check the length of the isoprenoid chain (7 to 25 isoprene units)
    chain_length = len(isoprenoid_matches[0])
    min_length = 7 * 5  # 7 isoprene units, each with 5 atoms
    max_length = 25 * 5  # 25 isoprene units, each with 5 atoms

    if chain_length < min_length or chain_length > max_length:
        return False, "Isoprenoid chain length outside the typical range for polyprenol phosphates"

    # Check if the isoprenoid chain is attached to the phosphate group via a terminal double bond
    for match in isoprenoid_matches:
        start_atom = mol.GetAtomWithIdx(match[0])
        end_atom = mol.GetAtomWithIdx(match[-1])

        if start_atom.IsInRing() or end_atom.IsInRing():
            continue  # Skip if part of a ring

        if end_atom.HasBondToAtom(phosphate_pattern):
            return True, "Contains a linear isoprenoid chain attached to a phosphate group"

    return False, "Isoprenoid chain not attached to phosphate group"