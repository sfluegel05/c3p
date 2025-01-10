"""
Classifies: CHEBI:61655 steroid saponin
"""
from rdkit import Chem

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin typically has a steroid backbone with hydroxyl groups and attached sugar moieties.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a steroid saponin, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Generalized SMARTS pattern for steroid backbone
    # Four fused rings: three six-membered followed by a five-membered ring
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3CC4=CC=CC=4C3CC2C1")
    if steroid_pattern is None:
        return None, None  # Safety check if pattern creation failed
    
    # Check for steroid backbone
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for at least one hydroxyl group on the steroid backbone
    hydroxy_on_steroid_pattern = Chem.MolFromSmarts("C[OH]")  # Simplified
    if not mol.HasSubstructMatch(hydroxy_on_steroid_pattern):
        return False, "Steroid backbone lacks hydroxyl groups"
    
    # Generalize pattern for sugar moieties (ring structures with oxygen)
    sugar_ring_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1")
    sugar_matches = mol.GetSubstructMatches(sugar_ring_pattern)
    if len(sugar_matches) == 0:
        return False, "No sugar moieties attached"
    
    return True, "Contains a hydroxysteroid backbone with attached sugar moieties consistent with a steroid saponin"