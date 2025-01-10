"""
Classifies: CHEBI:23437 cyanohydrin
"""
from rdkit import Chem

def is_cyanohydrin(smiles: str):
    """
    Determines if a molecule is a cyanohydrin based on its SMILES string.
    A cyanohydrin is defined as an alpha-hydroxynitrile resulting from the formal addition of hydrogen cyanide to the C=O bond of an aldehyde or ketone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyanohydrin, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for cyanohydrin: a carbon with hydroxyl group and nitrile (C#N) attached
    cyanohydrin_pattern = Chem.MolFromSmarts("[C](O)[CX1]#N")
    
    # Check if the molecule matches the cyanohydrin pattern
    if mol.HasSubstructMatch(cyanohydrin_pattern):
        return True, "Molecule contains the defining substructure of a cyanohydrin"
    else:
        return False, "No cyanohydrin substructure found in the molecule"