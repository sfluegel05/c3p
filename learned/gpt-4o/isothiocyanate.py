"""
Classifies: CHEBI:52221 isothiocyanate
"""
from rdkit import Chem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is an isothiocyanate based on its SMILES string.
    An isothiocyanate contains the functional group N=C=S, often as a terminal group.

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
    
    # Look for isothiocyanate group pattern N=C=S, making it more specific to common bonds
    isothiocyanate_pattern = Chem.MolFromSmarts("N=C=S")
    
    # Find matches for the isothiocyanate group in the molecule
    matches = mol.GetSubstructMatches(isothiocyanate_pattern)

    if matches:
        # Iterate over matches and confirm the pattern
        for match in matches:
            nx, cx, sx = match
            if mol.GetAtomWithIdx(sx).GetDegree() == 1:  # Ensures sulfur is terminal
                return True, "Contains terminal isothiocyanate group (N=C=S)"
        
        return False, "Isothiocyanate group not terminal"
    else:
        return False, "No isothiocyanate group (N=C=S) found"