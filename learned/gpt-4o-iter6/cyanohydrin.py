"""
Classifies: CHEBI:23437 cyanohydrin
"""
from rdkit import Chem

def is_cyanohydrin(smiles: str):
    """
    Determines if a molecule is a cyanohydrin based on its SMILES string.
    A cyanohydrin contains an alpha-hydroxynitrile group resulting from the addition of HCN
    to the C=O bond of an aldehyde or ketone.

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

    # Define the corrected cyanohydrin pattern
    # A carbon with a hydroxyl group (OH) and an adjacent nitrile group (C#N)
    cyanohydrin_pattern = Chem.MolFromSmarts("[#6]([#8H1])[#6]#[#7]")

    if mol.HasSubstructMatch(cyanohydrin_pattern):
        return True, "Contains a hydroxyl group adjacent to a nitrile group"
    else:
        return False, "Does not contain the cyanohydrin functional group"