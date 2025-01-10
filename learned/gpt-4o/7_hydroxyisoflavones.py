"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
from rdkit import Chem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    A 7-hydroxyisoflavone has a hydroxy group at the 7-position on the isoflavone skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 7-hydroxyisoflavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the isoflavone pattern, focusing on the core structure with a hydroxy group at position 7
    # The pattern should capture the base structure of 7-hydroxyisoflavones
    seven_hydroxyisoflavone_pattern = Chem.MolFromSmarts("c1cc2occ(c2c(=O)c3ccc(O)cc13)")

    # Check if the molecule has the 7-hydroxyisoflavone pattern
    if mol.HasSubstructMatch(seven_hydroxyisoflavone_pattern):
        return True, "Contains the 7-hydroxyisoflavone core structure"
    else:
        return False, "Does not have the 7-hydroxyisoflavone core structure"