"""
Classifies: CHEBI:48953 cyclohexenones
"""
from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is typically a cyclohexene ring with a ketone group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define cyclohexenone core pattern: six-membered ring with a ketone
    cyclohexenone_pattern = Chem.MolFromSmarts("C1=CC(=O)CCC1")
    
    # Check if cyclohexenone pattern is present
    if mol.HasSubstructMatch(cyclohexenone_pattern):
        return True, "Contains cyclohexenone core (six-membered ring with ketone group)"
    
    return False, "Cyclohexenone core not found"

# This function can be used to classify molecules based on the presence of a cyclohexenone structure.