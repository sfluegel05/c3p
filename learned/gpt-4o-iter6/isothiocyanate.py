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
    
    # Look for the isothiocyanate group with the R-group connection
    # Pattern updated to check that N is bonded to at least one carbon (R group)
    isothiocyanate_pattern = Chem.MolFromSmarts("[#6][#7]=[#6]=[#16]")  # Organic R connected to N=C=S
    matches = mol.GetSubstructMatches(isothiocyanate_pattern)
    if not matches:
        # Optional: fallback simple check (only if specific R-pattern check fails)
        isothiocyanate_fallback_pattern = Chem.MolFromSmarts("[#7]=[#6]=[#16]")
        matches = mol.GetSubstructMatches(isothiocyanate_fallback_pattern)

    # If matches are found, check the length/networks of carbon chain
    if matches:
        for match in matches:
            nitrogen = match[1]
            nitrogen_neighbours = mol.GetAtomWithIdx(nitrogen).GetNeighbors()
            for neighbor in nitrogen_neighbours:
                if neighbor.GetAtomicNum() == 6:  # Ensure Carbon-R group
                    # Check if Carbon is part of a larger chain (representative R group)
                    carbon_chain_size = len(list(neighbor.GetNeighbors()))
                    if carbon_chain_size > 1:  # Example: considering R chains of at least 2 atoms
                        return True, "Contains isothiocyanate group (R-N=C=S)"

    return False, "Does not contain isothiocyanate group (R-N=C=S)"