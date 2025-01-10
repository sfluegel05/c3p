"""
Classifies: CHEBI:30527 flavin
"""
from rdkit import Chem

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin based on its SMILES string.
    A flavin is a derivative of the dimethylisoalloxazine core with a substitution at the 10 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a flavin, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Define the dimethylisoalloxazine-like core structure
    dimethylisoalloxazine_core = Chem.MolFromSmarts("Cc1c(C)cc2nc3c([nH]c(=O)[nH]c3=O)n2cc1")
    
    # Check if the molecule contains the core structure
    core_matches = mol.GetSubstructMatches(dimethylisoalloxazine_core)
    if not core_matches:
        return False, "Does not contain the dimethylisoalloxazine core"
    
    # Analyze each match to find substitution at the correct position
    for match in core_matches:
        # Determine the position corresponding to '10' in the core
        core_at_index = match[5]  # Adjust index based on manual validation if 5 is '10'
        
        # Get the atom and its neighbors
        core_atom = mol.GetAtomWithIdx(core_at_index)
        
        # Check if there is a substitution (at least one neighbor not in core)
        is_substituted = any(neighbor.GetIdx() not in match for neighbor in core_atom.GetNeighbors())
        
        if is_substituted:
            return True, "Contains dimethylisoalloxazine core with substitution at position 10"

    return False, "No substitution at position 10 of the dimethylisoalloxazine core"