"""
Classifies: CHEBI:52221 isothiocyanate
"""
from rdkit import Chem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is an isothiocyanate based on its SMILES string.
    Isothiocyanates have the functional group N=C=S.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isothiocyanate, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to create an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Simple SMARTS pattern for isothiocyanate group: N=C=S
    isothiocyanate_pattern = Chem.MolFromSmarts("N=C=S")

    # Check for the presence of the isothiocyanate functional group
    if mol.HasSubstructMatch(isothiocyanate_pattern):
        return True, "Contains isothiocyanate functional group (N=C=S)"
    else:
        return False, "Does not contain isothiocyanate functional group (N=C=S)"