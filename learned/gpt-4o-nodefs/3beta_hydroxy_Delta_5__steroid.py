"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem

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

    # SMARTS pattern for the steroid core
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CC3C4C2CCC4CCC3")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Steroid backbone not found"

    # Check for 3beta-hydroxy group
    hydroxy_3beta_pattern = Chem.MolFromSmarts("[C@@H]1([C@H]2CCCC3CC(O)CCC3C2)C[C@@H](C1)O")
    if not mol.HasSubstructMatch(hydroxy_3beta_pattern):
        return False, "3beta-hydroxy group not found"

    # Check for 5-ene (double bond in Delta(5)-position)
    delta5_double_bond_pattern = Chem.MolFromSmarts("C1=C[C@@H]([C@@H]2CCCC3CCCCC3C2)C[C@@H](C1)O")
    if not mol.HasSubstructMatch(delta5_double_bond_pattern):
        return False, "Delta(5)-double bond not found"

    return True, "Molecule is classified as a 3beta-hydroxy-Delta(5)-steroid"