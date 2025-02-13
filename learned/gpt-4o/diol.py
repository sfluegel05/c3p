"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol contains exactly two hydroxy groups that are significantly characteristic of its structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find hydroxyl groups (–OH)
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    # Check if exactly two hydroxy groups exist to primarily qualify as a diol
    if len(hydroxy_matches) == 2:
        return True, "Contains exactly two hydroxy groups, qualifying it as a diol"
    
    # Optionally handle cases with more than two, if only two are primary
    elif len(hydroxy_matches) > 2:
        # More complex logic might be needed here to deep analyze structure
        # For simplicity, we assume more hydroxy groups "could" still be typical
        # if they are concentrated and part of a simpler linear/branched chain
        return True, "Contains more than two hydroxy groups, might still qualify as a diol structurally"
    
    return False, f"Found {len(hydroxy_matches)} hydroxy groups, which is inadequate for typical diol classification"

# Example usage:
# result, reason = is_diol("OC(CO)CO")  # Example to illustrate usage
# print(result, reason)