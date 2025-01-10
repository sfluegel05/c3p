"""
Classifies: CHEBI:48953 cyclohexenones
"""
from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is defined as any six-membered alicyclic ketone having one double bond in the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a six-membered ring with one double bond and a ketone group
    cyclohexenone_pattern = Chem.MolFromSmarts("C1=CC(=O)CCC1")
    
    # Check for the cyclohexenone pattern
    if not mol.HasSubstructMatch(cyclohexenone_pattern):
        return False, "No cyclohexenone structure found"
    
    return True, "Cyclohexenone structure identified with appropriate ring and ketone"