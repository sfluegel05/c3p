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
    dimethylisoalloxazine_core = Chem.MolFromSmarts("Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C)cc2cc1C")

    # Check if the molecule contains the core structure
    core_match = mol.GetSubstructMatch(dimethylisoalloxazine_core)
    if not core_match:
        return False, "Does not contain the dimethylisoalloxazine core"

    # Mapping atom indexes in core_match to check the 10th position
    # Atom indices should reflect the relevant positions based on the SMARTS pattern
    carbon_10_index = core_match[7]  # Adjust index based on true position of 10 in matched SMARTS

    # Evaluate if there is a substitution at position 10
    carbon_10_atom = mol.GetAtomWithIdx(carbon_10_index)
    is_substituted = any(neighbor.GetIdx() not in core_match for neighbor in carbon_10_atom.GetNeighbors())
    
    if is_substituted:
        return True, "Contains dimethylisoalloxazine core with substitution at position 10"

    return False, "No substitution at position 10 of the dimethylisoalloxazine core"