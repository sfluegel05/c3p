"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is defined as a 3-hydroxy steroid closely related to cholestan-3-ol.

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

    # Steroid backbone SMARTS pattern: four-ring steroid core 
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4CCCC(C4)C3C2C1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for hydroxyl group at the 3-position (3-hydroxy group)
    steroid_3_oh_pattern = Chem.MolFromSmarts("C1(CCC2C3CCC4CCCC(C4)C3C2C1)O")
    if not mol.HasSubstructMatch(steroid_3_oh_pattern):
        return False, "No 3-hydroxy group found in the steroid backbone"

    return True, "Contains steroid backbone with 3-hydroxy group"