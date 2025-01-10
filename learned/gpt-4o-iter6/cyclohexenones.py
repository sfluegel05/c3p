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
    
    # Define a SMARTS pattern for a six-membered ring with one double bond
    cyclohexene_pattern = Chem.MolFromSmarts("C1=CCCCC1")
    
    # Check for the six-membered ring with one double bond
    if not mol.HasSubstructMatch(cyclohexene_pattern):
        return False, "No six-membered cyclohexene ring found"
    
    # Define a SMARTS pattern for a ketone group within a ring
    ketone_in_ring_pattern = Chem.MolFromSmarts("C1(=O)[#6]1")
    
    # Check for a ketone group in the ring
    if not mol.HasSubstructMatch(ketone_in_ring_pattern):
        return False, "No ketone group in the cyclohexene ring found"
    
    return True, "Cyclohexenone structure identified with a cyclohexene ring containing a ketone group"