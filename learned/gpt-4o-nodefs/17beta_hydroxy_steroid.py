"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
from rdkit import Chem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    A 17beta-hydroxy steroid contains a steroid backbone with a hydroxyl group at the 17-beta position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a generalized SMARTS pattern for the steroid backbone (rings A, B, C, and D)
    steroid_skeleton_pattern = Chem.MolFromSmarts("C1CC[C@H]2C(C1)CCC3CCC4C3CCC4C2")
    if not mol.HasSubstructMatch(steroid_skeleton_pattern):
        return False, "No steroid skeleton found"

    # Define SMARTS pattern for hydroxyl group at the 17-beta position (more broadly)
    hydroxyl_17beta_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3(CCC(O)C4)")
    if not mol.HasSubstructMatch(hydroxyl_17beta_pattern):
        return False, "No hydroxyl group at 17-beta position found"

    return True, "Molecule matches 17beta-hydroxy steroid structure"