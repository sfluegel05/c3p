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
        bool: True if the molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for cyclohexenone structure
    # Matches six-membered ring with exactly one double bond and one ketone
    cyclohexenone_pattern = Chem.MolFromSmarts("C1=CC(=O)CCC1")
    
    # Check for the cyclohexenone pattern
    if mol.HasSubstructMatch(cyclohexenone_pattern):
        return True, "Cyclohexenone structure identified in the molecule"

    return False, "No cyclohexenone structure found"