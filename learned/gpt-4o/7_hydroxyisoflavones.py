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

    # Correct pattern for isoflavone core with 7-hydroxy group (ensures correct structural alignment)
    isoflavone_core_7_oh_pattern = Chem.MolFromSmarts("Oc1ccc2c(c1)c(=O)cc3c2occ3-c1ccccc1")

    # Check for the isoflavone core structure with a 7-hydroxy group
    if not mol.HasSubstructMatch(isoflavone_core_7_oh_pattern):
        return False, "No 7-hydroxyisoflavone structure found"
    
    return True, "Molecule is a 7-hydroxyisoflavone with the correct structure"