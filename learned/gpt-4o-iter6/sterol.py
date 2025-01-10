"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is defined as a 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for more flexible steroid-like backbone:
    # A pattern for a saturated polycyclic backbone with 3 six-membered rings and 1 five-membered ring
    # Allow for flexibility with unsaturation (double bonds)
    steroid_pattern = Chem.MolFromSmarts("C1(C)C[C@H]2[C@@H]3C=CC[C@H]4C(C)=C/C=C/[C@]4(C)C3CCC2C1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No flexible steroid-like backbone found"

    # Check for a hydroxyl group [-OH] which is characteristic for sterols
    # Consider hydroxyl groups attached to both aromatic and aliphatic carbons
    hydroxyl_pattern = Chem.MolFromSmarts("[#6][OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # Check for the presence of typical modifications or side chains present in sterols
    # Using a lenient approach here to account for variety
    possible_sidechain_pattern = Chem.MolFromSmarts("[CX3,CX4]C([CX3,CX4])")
    if mol.HasSubstructMatch(possible_sidechain_pattern):
        return True, "Contains steroid-like backbone with hydroxyl group and compatible side chains, consistent with sterol definition"
    else:
        # Even if side chain is not found, it might still be a sterol
        return True, "Contains steroid-like backbone with hydroxyl group; possible sterol with unconventional side chain"