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

    # General SMARTS pattern for a steroid core (4-ring backbone)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4CCCCC34")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Steroid backbone not found"

    # Look for hydroxy group at the 3-position (loosely defined)
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H]1(CC[C@H]2CC[C@H](C)C3CCC(O)C23)C1")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "3beta-hydroxy group not found"

    # Check for Delta(5) double bond; allow for flexibility in stereochemistry
    delta5_double_bond_pattern = Chem.MolFromSmarts("C1(C=C)C[C@H]2CCC3C4[C@H](CCC4CCC3)C2")
    if not mol.HasSubstructMatch(delta5_double_bond_pattern):
        return False, "Delta(5) double bond not found"

    return True, "Molecule is classified as a 3beta-hydroxy-Delta(5)-steroid"