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

    # SMARTS for steroid core (4 fused rings, at least two of them 6-membered)
    steroid_core_pattern = Chem.MolFromSmarts("[R6]1[R6]2[R]3[R]1[R]4[R]2[R]3[R6]4")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Molecule does not have a steroid core"

    # SMARTS for carbonyl at position 3. Carbonyl group should be directly attached to a carbon of a 6 membered ring
    oxo_at_3_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[C;R6]")
    if not mol.HasSubstructMatch(oxo_at_3_pattern):
        return False, "Molecule does not have a carbonyl at position 3"


    return True, "Molecule is a 3-oxo steroid"