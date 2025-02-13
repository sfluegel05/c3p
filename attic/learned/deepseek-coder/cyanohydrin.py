"""
Classifies: CHEBI:23437 cyanohydrin
"""
"""
Classifies: CHEBI:51853 cyanohydrin
"""
from rdkit import Chem

def is_cyanohydrin(smiles: str):
    """
    Determines if a molecule is a cyanohydrin based on its SMILES string.
    A cyanohydrin is an alpha-hydroxynitrile resulting from the formal addition of hydrogen cyanide to the C=O bond of an aldehyde or ketone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyanohydrin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the cyanohydrin pattern: a carbon with both -OH and -C#N groups
    cyanohydrin_pattern = Chem.MolFromSmarts("[CX4]([OH])[C]#N")
    if not mol.HasSubstructMatch(cyanohydrin_pattern):
        return False, "No cyanohydrin pattern found (carbon with both -OH and -C#N groups)"

    return True, "Contains a carbon with both -OH and -C#N groups, characteristic of a cyanohydrin"