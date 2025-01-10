"""
Classifies: CHEBI:68452 azole
"""
from rdkit import Chem

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole based on its SMILES string.
    Azoles are defined as a monocyclic heteroarene consisting of a five-membered ring containing nitrogen.
    The ring may also contain one or more other non-carbon atoms such as nitrogen, sulfur, or oxygen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an azole, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for azole ring: five-membered aromatic ring with at least one nitrogen atom
    # Additional non-carbon atoms allowed are O or S
    azole_pattern = Chem.MolFromSmarts("n1[c,o,s,n;R]c[c,o,s,n;R]c1")

    # Check if the molecule has a substructure match with the azole pattern
    if mol.HasSubstructMatch(azole_pattern):
        return True, "Contains a five-membered heteroarene ring with nitrogen"

    return False, "Does not contain a five-membered heteroarene ring with nitrogen"