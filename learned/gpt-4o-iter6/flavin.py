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

    # Define a more comprehensive dimethylisoalloxazine core pattern
    dimethylisoalloxazine_core = Chem.MolFromSmarts("CC1=CC2=C(C=C1C)N3C(=O)NC(=O)C3=NC2") 
    # Mapping adjustment for expected different substitutions
    core_matches = mol.GetSubstructMatches(dimethylisoalloxazine_core)
    
    if not core_matches:
        return False, "Does not contain the dimethylisoalloxazine core"

    for match in core_matches:
        # Assuming the 10 position is roughly equivalent to the 5th index in this pattern
        pos10_atom_idx = match[5]
        pos10_atom = mol.GetAtomWithIdx(pos10_atom_idx)
        
        # Check if there is a substitution at position '10'
        if len(pos10_atom.GetNeighbors()) > 1:  # More than core's direct bonding
            return True, "Contains dimethylisoalloxazine core with substitution at position 10"

    return False, "No substitution at position 10 of the dimethylisoalloxazine core"