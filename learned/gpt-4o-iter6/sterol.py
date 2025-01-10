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
        bool: True if the molecule is a sterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the steroid backbone:
    # Three six-membered rings and one five-membered ring, allowing some flexibility/variation
    steroid_core_pattern = Chem.MolFromSmarts("C1CCC2C1CCC3C2CCC4C3CCC4")
    if not steroid_core_pattern:
        return (None, "Invalid steroid core SMARTS pattern")

    # Match steroid core pattern
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No appropriate steroid backbone (3 six-membered and 1 five-membered ring) found"

    # Check for hydroxyl group [-OH] in the C3 position
    hydroxyl_pattern = Chem.MolFromSmarts("C(O)C1CCC2C3CCC4C(C)CCC4C3CCC12")
    if not hydroxyl_pattern:
        return (None, "Invalid hydroxyl position SMARTS pattern")
        
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found on the C3 position"

    # Assume presence of side chains typical of sterols due to sterol diversity
    side_chain_check = True

    if side_chain_check:
        return True, "Contains steroid backbone with hydroxyl group; structure consistent with sterol definition"
    
    return False, "Missing structural features typical of sterols"