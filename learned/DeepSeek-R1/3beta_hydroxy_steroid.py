"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: 3beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    A 3beta-hydroxy steroid has a hydroxyl group in the beta position at carbon 3 of the steroid nucleus.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for steroid nucleus using a SMARTS pattern for the tetracyclic structure
    steroid_smarts = Chem.MolFromSmarts("[C]1([C]2([C]3([C]([C]4([C]([C]([C]([C]1)[C]([C]2)[C]3)[C]4)))))")
    if not mol.HasSubstructMatch(steroid_smarts):
        return False, "No steroid nucleus found"

    # Check for hydroxyl group at position 3 with beta configuration
    # Using SMARTS to find -OH on carbon 3 with beta (R) configuration
    beta_oh_smarts = Chem.MolFromSmarts("[C@H]1([OH])CC[C@H]2[C@@H]3CCC4CC[C@H](C4)[C@H]3CC[C@]12C")
    if mol.HasSubstructMatch(beta_oh_smarts):
        return True, "3beta-hydroxy group present on steroid nucleus"

    return False, "No 3beta-hydroxy group found"