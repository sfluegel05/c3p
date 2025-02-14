"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    A 16beta-hydroxy steroid has a hydroxyl group at the 16th carbon of a steroid core,
    where this hydroxy group has a beta-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 16beta-hydroxy steroid, False otherwise
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a flexible steroid core SMARTS pattern with a beta-hydroxy at C16
    # The pattern specifies the 4 fused rings (3x 6-membered, 1x 5-membered) and
    # specifically looks for a hydroxyl group with a beta configuration at C16
    # [C] is any sp3 carbon, and [C@H] specifies that the carbon atom
    # has a specific stereochemistry (clockwise), which defines it as beta.
    steroid_smarts_beta_oh = "[C]1[C][C]2[C]3[C]([C]1)[C][C]4[C]3[C]([C]2)[C]([C]5)[C]4([C])[C]([C@H]5O)"
    steroid_pattern_beta_oh = Chem.MolFromSmarts(steroid_smarts_beta_oh)

    # Find substructure match
    match = mol.HasSubstructMatch(steroid_pattern_beta_oh)
    if not match:
         return False, "No 16-beta hydroxyl steroid core found"

    return True, "16beta-hydroxy steroid found"