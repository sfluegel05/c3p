"""
Classifies: CHEBI:48953 cyclohexenones
"""
from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is a six-membered alicyclic ketone with one double bond in the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string into a molecular representation
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create a more generic SMARTS pattern for a six-membered ring with one double bond and a ketone
    cyclohexenone_pattern = Chem.MolFromSmarts("C1=CC(=O)CCC1 |$;c;;;;;c$|")
    
    # Check if the molecule matches the pattern
    if not mol.HasSubstructMatch(cyclohexenone_pattern):
        return False, "Does not contain a six-membered ring with one double bond and a ketone"

    return True, "Contains a six-membered alicyclic ketone with one double bond in the ring"