"""
Classifies: CHEBI:52221 isothiocyanate
"""
"""
Classifies: Isothiocyanate
"""
from rdkit import Chem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is an isothiocyanate based on its SMILES string.
    An isothiocyanate is characterized by the -N=C=S group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isothiocyanate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for isothiocyanate group
    isothiocyanate_pattern = Chem.MolFromSmarts("[N]=[C]=[S]")

    # Check for substructure match
    if mol.HasSubstructMatch(isothiocyanate_pattern):
        return True, "Contains isothiocyanate group"
    else:
        return False, "Missing isothiocyanate group"