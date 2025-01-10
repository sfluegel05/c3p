"""
Classifies: CHEBI:16385 organic sulfide
"""
from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide based on its SMILES string.
    An organic sulfide contains a sulfur atom bonded to two organic groups (R-S-R' with R =/= H).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an organic sulfide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # sulfur atom that is singly-bonded to two carbon atoms
    sulfide_pattern = Chem.MolFromSmarts("C-S-C")
    if mol.HasSubstructMatch(sulfide_pattern):
        return True, "Contains the R-S-R' structure characteristic of organic sulfides"
    
    return False, "Does not contain the R-S-R' structure typically found in organic sulfides"