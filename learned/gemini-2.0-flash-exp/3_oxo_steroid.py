"""
Classifies: CHEBI:47788 3-oxo steroid
"""
from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is a steroid with a carbonyl group (=O) at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for the steroid core with specific ring fusion
    # Based on the generic steroid structure:
    #    1
    #   / \
    #  2---3
    # / \ / \
    # 4---5---6
    # |   |   |
    # 7---8---9
    #  \ / \ /
    #   10--11
    #     \ /
    #      12
    # Ring A: 1-2-4-7-10-11
    # Ring B: 2-3-5-8-10-11
    # Ring C: 3-5-6-9-8
    # Ring D: 6-9-12-11
    steroid_core_pattern = Chem.MolFromSmarts("[C;R6]1[C;R6]2[C;R6]3[C;R5]4[C;R6]1[C;R6]2[C;R6]3[C;R5]4")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Molecule does not have a steroid core"

    # SMARTS for carbonyl at position 3 (specifically, at the second carbon of ring A)
    oxo_at_3_pattern = Chem.MolFromSmarts("[C;R6]1[C;R6](=[OX1])[C;R6]2[C;R5]3[C;R6]1[C;R6]2[C;R6]3")
    if not mol.HasSubstructMatch(oxo_at_3_pattern):
        return False, "Molecule does not have a carbonyl at position 3"

    return True, "Molecule is a 3-oxo steroid"