"""
Classifies: CHEBI:74716 withanolide
"""
"""
Classifies: CHEBI:51800 withanolide

A withanolide is any steroid lactone that is a C28 steroid with a modified side chain
forming a lactone ring and its substituted derivatives.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_withanolide(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a withanolide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid backbone
    steroid_pattern = Chem.MolFromSmarts("[C@]1(C)(CC[C@]2([H])[C@@]3([H])CC[C@@]4([H])C[C@H](C)[C@@H](O)[C@@H](C)[C@H](O)[C@@H](C)[C@@H](O)[C@@H](C)[C@@H]4CCC3=C2)CC1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Does not contain steroid backbone"

    # Check for modified side chain with lactone ring
    sidechain_pattern = Chem.MolFromSmarts("[C@@]1([C@@](CC1=O)(C)[C@@]2([H])CC[C@@]3([H])[C@]4([H])CC[C@@]5([H])C[C@H](C)[C@@H](O)[C@@H](C)[C@H](O)[C@@H](C)[C@@H]5CCC4=C3C2)C")
    if not mol.HasSubstructMatch(sidechain_pattern):
        return False, "No modified side chain with lactone ring found"

    # Count number of rings (withanolides typically have 4-5 rings)
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()
    if n_rings < 4 or n_rings > 5:
        return False, "Number of rings outside typical range for withanolides"

    return True, "Contains steroid backbone with modified side chain forming a lactone ring"