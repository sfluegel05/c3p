"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: methyl sulfide
"""

from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is any molecule where a sulfur atom is connected via single bonds
    to at least one methyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for methyl sulfide
    methyl_sulfide_pattern = Chem.MolFromSmarts("[#16X2]([CH3])[*]")
    
    # Search for methyl sulfide pattern
    if mol.HasSubstructMatch(methyl_sulfide_pattern):
        return True, "Contains methyl sulfide group"
    else:
        return False, "Does not contain methyl sulfide group"