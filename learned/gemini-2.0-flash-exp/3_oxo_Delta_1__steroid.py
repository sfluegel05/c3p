"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1)-steroid based on its SMILES string.
    A 3-oxo-Delta(1)-steroid is a steroid with a ketone at position 3 and a double bond
    between positions 1 and 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(1)-steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the basic steroid core structure - 4 rings
    # The 'A' allows for any atom but H
    steroid_core_pattern = Chem.MolFromSmarts("[A]1[A]2[A]3[A]4[A]1[A]5[A]2[A]3[A]45") # general steroid core pattern.

    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Molecule does not have a steroid core."

   # Define a SMARTS pattern for the 3-oxo group. Note this must be attached to a ring.
    oxo_3_pattern = Chem.MolFromSmarts("[C]1(=O)[C]([A])[A]2[A]3[A]4[A]1[A]5[A]2[A]3[A]45")

    if not mol.HasSubstructMatch(oxo_3_pattern):
        return False, "Molecule does not have a 3-oxo group"

    # Define a SMARTS pattern for the double bond between C1 and C2
    delta_1_pattern = Chem.MolFromSmarts("[C]1=[C][C]2[A]3[A]4[A]1[A]5[A]2[A]3[A]45") # double bond from position 1 to position 2.

    if not mol.HasSubstructMatch(delta_1_pattern):
        return False, "Molecule does not have a double bond between positions 1 and 2."
    
    
    return True, "Molecule is a 3-oxo-Delta(1)-steroid."