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
    
    # Look for the specific isothiocyanate functional group with appropriate connectivity
    isothiocyanate_pattern = Chem.MolFromSmarts("[#7]=[#6]=[#16]")  # Specific pattern for N=C=S
    matches = mol.GetSubstructMatches(isothiocyanate_pattern)
    if matches:
        # Additional check: ensure N is connected to a carbon chain (R group)
        for match in matches:
            nitrogen = match[0]
            nitrogen_neighbours = mol.GetAtomWithIdx(nitrogen).GetNeighbors()
            for neighbor in nitrogen_neighbours:
                if neighbor.GetAtomicNum() == 6:  # Looking for carbon attached to Nitrogen
                    return True, "Contains isothiocyanate group (R-N=C=S)"
    
    return False, "Does not contain isothiocyanate group (R-N=C=S)"