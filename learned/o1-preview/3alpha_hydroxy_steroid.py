"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    A 3alpha-hydroxy steroid is a steroid with a hydroxyl group at position 3 in the alpha orientation.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid nucleus SMARTS pattern (cyclopentanoperhydrophenanthrene)
    steroid_smarts = """
    [#6]1([#6])[#6][#6]2[#6]([#6]1)[#6][#6]3[#6]([#6]2)[#6][#6]4[#6]([#6]3)[#6][#6][#6][#6]4
    """
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Steroid nucleus not found"

    # Define 3alpha-hydroxy group pattern on steroid nucleus
    # The numbering corresponds to the steroid backbone
    # We need to identify the 3-position carbon with alpha-oriented hydroxyl group
    # Since RDKit does not inherently provide atom numbering, we map the pattern
    # to the molecule and check the stereochemistry

    # SMARTS pattern for 3alpha-hydroxy steroid
    alpha_hydroxyl_smarts = """
    [#6]1([#6])[#6]([C@H](O))[#6]2[#6]([#6]1)[#6][#6]3[#6]([#6]2)[#6][#6]4[#6]([#6]3)[#6][#6][#6][#6]4
    """
    alpha_hydroxyl_pattern = Chem.MolFromSmarts(alpha_hydroxyl_smarts)
    matches = mol.GetSubstructMatches(alpha_hydroxyl_pattern, useChirality=True)

    if matches:
        return True, "Contains 3alpha-hydroxy group with correct stereochemistry"
    else:
        # Check if hydroxyl group is present at position 3 but with incorrect stereochemistry
        beta_hydroxyl_smarts = """
        [#6]1([#6])[#6]([C@@H](O))[#6]2[#6]([#6]1)[#6][#6]3[#6]([#6]2)[#6][#6]4[#6]([#6]3)[#6][#6][#6][#6]4
        """
        beta_hydroxyl_pattern = Chem.MolFromSmarts(beta_hydroxyl_smarts)
        beta_matches = mol.GetSubstructMatches(beta_hydroxyl_pattern, useChirality=True)
        if beta_matches:
            return False, "Hydroxyl group at position 3 has beta orientation (not alpha)"
        else:
            return False, "No 3alpha-hydroxy group found"