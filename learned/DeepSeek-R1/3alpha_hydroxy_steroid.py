"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
"""
Classifies: 3alpha-hydroxy steroid (CHEBI:36807)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    A 3alpha-hydroxy steroid has a hydroxyl group in the alpha position at carbon 3 of the steroid nucleus.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for steroid nucleus (four fused rings: three 6-membered, one 5-membered)
    # Simplified check using a SMARTS pattern for the steroid core
    steroid_core = Chem.MolFromSmarts("[C]1[C@@H]2[C@H](C3[C@@H](C[C@H]4CCCC34)CC2)CC1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid nucleus detected"

    # Check for hydroxyl group at position 3 with alpha configuration
    # SMARTS pattern for 3alpha-OH: hydroxyl attached to C3 with specific stereochemistry
    alpha_oh_pattern = Chem.MolFromSmarts("[C@@H]3([C@H](CC1)CC2)[OH]")
    if not mol.HasSubstructMatch(alpha_oh_pattern):
        return False, "No 3alpha-hydroxy group found"

    return True, "3alpha-hydroxy group present on steroid nucleus"