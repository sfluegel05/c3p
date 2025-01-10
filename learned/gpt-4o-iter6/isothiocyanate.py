"""
Classifies: CHEBI:52221 isothiocyanate
"""
from rdkit import Chem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is an isothiocyanate based on its SMILES string.
    An isothiocyanate contains a functional group with the general formula R-N=C=S.

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

    # Look for the standalone isothiocyanate pattern (R-N=C=S)
    # We use a more generic pattern adjusted to accommodate more linkage varieties
    isothiocyanate_pattern = Chem.MolFromSmarts("[$([NX2]=C=S)]")  # Generic pattern for N=C=S
    if mol.HasSubstructMatch(isothiocyanate_pattern):
        return True, "Contains isothiocyanate group (R-N=C=S)"

    return False, "Does not contain isothiocyanate group (R-N=C=S)"