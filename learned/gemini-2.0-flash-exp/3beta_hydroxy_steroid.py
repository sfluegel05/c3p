"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: 3beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    A 3beta-hydroxy steroid is a steroid with a hydroxyl group in the beta position on the 3rd carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a simplified SMARTS pattern to capture the essential steroid core and the beta-hydroxyl.
    # The pattern below defines a four-ring core with a beta hydroxyl at the 3rd position, and uses a SMARTS atom map for verification.
    # [C@H] specifies beta configuration.
    # [O;H1] will verify that this is a hydroxyl group with one implied H atom.
    steroid_pattern = Chem.MolFromSmarts("[C]12[C]3[C]4[C]([C]([C]1[C]2)CC3)CC[C@H]([O;H1])4")

    if steroid_pattern is None:
      return False, "Invalid SMARTS pattern"

    # Check if the molecule matches the pattern
    matches = mol.GetSubstructMatches(steroid_pattern)

    if not matches:
        return False, "Molecule does not match the steroid core with beta-hydroxyl at position 3"
    
    return True, "Molecule matches the criteria for a 3beta-hydroxy steroid"