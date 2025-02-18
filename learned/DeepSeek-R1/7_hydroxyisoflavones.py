"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
"""
Classifies: CHEBI:7-hydroxyisoflavones
"""
from rdkit import Chem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    A 7-hydroxyisoflavone has a hydroxy group at the 7-position of the isoflavone core structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 7-hydroxyisoflavone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # SMARTS pattern for 7-hydroxyisoflavone core:
    # Hydroxyl attached to benzene ring (position 7) fused to pyrone with oxygen-linked substituent
    pattern = Chem.MolFromSmarts("[OH]-c1ccc2c(c1)oc(-*)cc2=O")
    if mol.HasSubstructMatch(pattern):
        return True, "7-hydroxy group present on isoflavone core"
    
    # Alternative check for different substituents on the B ring
    # Broader pattern allowing any substitution on the B ring
    # Check if any hydroxyl is attached to the A ring at position 7
    # This approach may require more complex atom mapping
    
    return False, "No 7-hydroxy group detected on isoflavone core"