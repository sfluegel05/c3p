"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy-Delta(5)-steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the steroid core (4 fused rings, 3 six membered and 1 five membered)
    steroid_core_pattern = Chem.MolFromSmarts("[r6]12[r6][r5][r6]1[r6][r6][r6]2")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Steroid core not found"

    # Define SMARTS pattern for the 3-beta-hydroxyl group.
    beta_hydroxy_pattern = Chem.MolFromSmarts("[C@H]([OH])[C]([C])")
    if not mol.HasSubstructMatch(beta_hydroxy_pattern):
        return False, "No 3-beta hydroxyl group found"

    # Define SMARTS pattern for the double bond between C5 and C6. The carbons 5 and 6 must belong to the steriod skeleton
    delta5_bond_pattern = Chem.MolFromSmarts("[C]1[C](=[C])[C]([C]2[C]1[C]3[C]2[C])3")
    if not mol.HasSubstructMatch(delta5_bond_pattern):
         return False, "No double bond between C5 and C6"

    return True, "3beta-hydroxy-Delta(5)-steroid identified"