"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    A 17beta-hydroxy steroid has the basic steroid core and a beta-configured hydroxyl group at position 17.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for the steroid core, improved to be less strict
    steroid_core_pattern = Chem.MolFromSmarts("[C:1]1[C:2][C:3][C:4]2[C:5][C:6]([C:1])([C:7][C:8]([C:2])C)C[C:3][C:4]C")

    match = mol.GetSubstructMatch(steroid_core_pattern)
    if not match:
        return False, "Molecule does not have the basic steroid core structure."

    # Get the C17 atom
    c17_idx = match[7] #Index 7 in the match array corresponds to the C17 atom

    # SMARTS pattern for beta-hydroxyl at C17 (explicit H)
    beta_oh_pattern_with_H = Chem.MolFromSmarts("[C@]([C])([C])[OH]")

    # SMARTS pattern for beta-hydroxyl at C17 (implicit H)
    beta_oh_pattern_no_H = Chem.MolFromSmarts("[C@]([C])([C])O")

    # Check for match to either of the patterns
    match_with_H = mol.GetSubstructMatch(beta_oh_pattern_with_H)
    match_no_H = mol.GetSubstructMatch(beta_oh_pattern_no_H)
    
    c17_match = None
    if match_with_H:
        for idx in match_with_H:
            if idx == c17_idx:
              c17_match = True
    if match_no_H:
        for idx in match_no_H:
           if idx == c17_idx:
              c17_match = True

    if not c17_match:
        return False, "Hydroxyl group at position 17 is not beta-configured"
    
    return True, "Molecule contains a steroid core with a beta-hydroxy group at position 17."