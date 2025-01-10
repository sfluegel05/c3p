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
    
    # Define a more accurate SMARTS pattern for the steroid scaffold (three 6-membered and one 5-membered rings fused)
    steroid_skeleton_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4(C3CCC4)C")
    if not mol.HasSubstructMatch(steroid_skeleton_pattern):
        return False, "No steroid skeleton found"
    
    # Define a SMARTS pattern for hydroxyl group at 17-beta position
    # 17 position would generally be flexible as for its placement relative to ring A.
    # [C@] and [C@H] indicate stereochemistry and hydrogen
    hydroxyl_17beta_pattern = Chem.MolFromSmarts("C[C@]12CC[C@H](C)[C@@H](O)[C@]1(CC=C2)C")
    if not mol.HasSubstructMatch(hydroxyl_17beta_pattern):
        return False, "No hydroxyl group at 17-beta position found"
    
    return True, "Molecule matches 17beta-hydroxy steroid structure"