"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: 3-oxo steroid
"""

from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is a steroid with a ketone (=O) functional group at position 3 of the steroid nucleus.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for steroid nucleus (cyclopentanoperhydrophenanthrene)
    steroid_nucleus_smarts = '[#6]12[#6][#6][#6]3([#6]1[#6][#6][#6]([#6]2)[#6]4[#6][#6][#6][#6][#6]34)'
    steroid_nucleus = Chem.MolFromSmarts(steroid_nucleus_smarts)
    if steroid_nucleus is None:
        return False, "Error in steroid nucleus SMARTS pattern"

    # Check if the molecule has the steroid nucleus
    if not mol.HasSubstructMatch(steroid_nucleus):
        return False, "Molecule does not contain steroid nucleus"

    # Define SMARTS pattern for ketone at position 3 of steroid nucleus
    # This pattern looks for a ketone (=O) attached to the third carbon in the A ring
    three_oxo_steroid_smarts = '''
    [#6]12[#6](=O)[#6][#6]3([#6]1[#6][#6][#6]([#6]2)[#6]4[#6][#6][#6][#6][#6]34)
    '''
    three_oxo_steroid_pattern = Chem.MolFromSmarts(three_oxo_steroid_smarts)
    if three_oxo_steroid_pattern is None:
        return False, "Error in 3-oxo steroid SMARTS pattern"

    # Check if molecule matches the 3-oxo steroid pattern
    if mol.HasSubstructMatch(three_oxo_steroid_pattern):
        return True, "Molecule is a 3-oxo steroid"
    else:
        return False, "Molecule does not have ketone at position 3 of steroid nucleus"