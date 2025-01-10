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
    
    # Define the isoflavone pattern with a hydroxy group specifically at the 7-position
    # The pattern covers flexibility for substitutions on the phenyl and heterocyclic rings without losing the core 7-hydroxy group feature
    isoflavone_core_pattern = Chem.MolFromSmarts("c1ccc2c(c1)oc(=O)c1c(ccc(O)c1)c2")
    
    # Check if the molecule has the 7-hydroxyisoflavone pattern
    if mol.HasSubstructMatch(isoflavone_core_pattern):
        return True, "Matches the core structure with a hydroxy group at the 7-position on the isoflavone framework"
    else:
        return False, "Does not match the core structure with a hydroxy group at the 7-position"