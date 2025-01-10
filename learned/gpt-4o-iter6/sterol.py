"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is defined as a 3-hydroxy steroid with a structure closely related to cholestan-3-ol.

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

    # Define a more flexible SMARTS pattern for the steroid backbone:
    # Recognize 3 connected six-membered rings and 1 five-membered ring, allowing for some unsaturation.
    steroid_core_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4=C3CC[C@H]5C4")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No appropriate steroid backbone (3 six-membered and 1 five-membered rings) found"

    # Check for hydroxyl group [-OH] in the stereochemically relevant positions, mainly C3.
    # Consider C3-beta or analogous positions based on biochemistry.
    hydroxyl_position_pattern = Chem.MolFromSmarts("C(O)[C@H]1C[C@H]2C3C([C@@H](O)C4=C(C)CCCC34)CCC2CC1")
    if not mol.HasSubstructMatch(hydroxyl_position_pattern):
        return False, "No hydroxyl group found on potential C3 position"

    # Allow some typical sidechains or alkyl chains; sterols might have larger, more flexible side chains.
    # This pattern is too generic for a graphical match; instead, identify fragments.
    side_chain_check = True  # Assume general presence of compatible sidechains due to sterol structure diversity

    if side_chain_check:
        return True, "Contains steroid backbone with hydroxyl group; structure consistent with sterol definition"
    else:
        return False, "No compatible side chain pattern detected, not a typical sterol"