"""
Classifies: CHEBI:23437 cyanohydrin
"""
from rdkit import Chem

def is_cyanohydrin(smiles: str):
    """
    Determines if a molecule is a cyanohydrin based on its SMILES string.
    A cyanohydrin contains a carbon atom attached to both a hydroxyl (-OH) and a nitrile (-C#N) group.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a cyanohydrin, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define cyanohydrin pattern as a carbon connected to both -OH and a -C#N group
    cyanohydrin_pattern = Chem.MolFromSmarts("[C]([OH])C#N")
    
    # Check for presence of cyanohydrin structural pattern
    if mol.HasSubstructMatch(cyanohydrin_pattern):
        return True, "Contains cyanohydrin moiety (C-OH and C#N)"
    else:
        return False, "Does not contain cyanohydrin moiety"

# Example usage:
# print(is_cyanohydrin("O[C@@H](C#N)c1ccccc1"))  # Should return (True, "...Contains cyanohydrin moiety...")