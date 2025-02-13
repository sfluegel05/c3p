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
    # Considering context or constraints around connections for specificity
    isothiocyanate_pattern = Chem.MolFromSmarts("[NX2]=[CX2]=[SX1]")  # basic pattern for N=C=S

    # A more strict pattern to ensure N=C=S is terminal or more isolated could be introduced
    # Additional constraints could be added here
    if mol.HasSubstructMatch(isothiocyanate_pattern):
        NCS_matches = mol.GetSubstructMatch(isothiocyanate_pattern)

        # Verify that the match is reasonable (e.g., hanging or certain connectivity)
        for match in NCS_matches:
            cx, nx, sx = match
            if mol.GetAtomWithIdx(sx).GetDegree() == 1:
                return True, "Contains terminal isothiocyanate group (N=C=S)"
        
        return False, "Isothiocyanate pattern present but not terminal or unique in structure"
    else:
        return False, "No isothiocyanate group (N=C=S) found"