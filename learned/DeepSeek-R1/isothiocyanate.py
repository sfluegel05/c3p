"""
Classifies: CHEBI:52221 isothiocyanate
"""
"""
Classifies: CHEBI:52221 isothiocyanate
"""
from rdkit import Chem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is an isothiocyanate based on its SMILES string.
    An isothiocyanate has the general formula R-N=C=S where R is any organic substituent.

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

    # Define SMARTS pattern for the isothiocyanate group: R-N=C=S
    # Matches any carbon attached to N which is double bonded to C=S
    isothiocyanate_pattern = Chem.MolFromSmarts("[C]-[N]=[C]=[S]")
    
    # Check for presence of the pattern
    if mol.HasSubstructMatch(isothiocyanate_pattern):
        return True, "Contains R-N=C=S group"
    else:
        return False, "No R-N=C=S group detected"