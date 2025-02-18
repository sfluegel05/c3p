"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
"""
Classifies: 7-hydroxyisoflavones
"""
from rdkit import Chem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    A 7-hydroxyisoflavone has a specific isoflavone structure with a hydroxyl group at the 7-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 7-hydroxyisoflavone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Isoflavone core pattern with a more flexible SMARTS
    # Consider variations that still conform to the basic structure
    isoflavone_core = Chem.MolFromSmarts("Oc1cc(-c2coCc3c2cc(O)cc3)cc(O)c1")  # Allow variations around the core
    
    # Check for the isoflavone core structure with hydroxyls, especially at position 7
    if not mol.HasSubstructMatch(isoflavone_core):
        return False, "No central 7-hydroxyisoflavone structure found"
    
    # Verify the exact position of hydroxyl groups can vary around the core
    hydroxyl_position_7 = Chem.MolFromSmarts("Oc1cc(-c2coCc3c2cc(O)cc3)ccccc1")
    if not mol.HasSubstructMatch(hydroxyl_position_7):
        return False, "No hydroxyl group found at the 7-position"

    return True, "Molecule is a 7-hydroxyisoflavone with the correct structure"