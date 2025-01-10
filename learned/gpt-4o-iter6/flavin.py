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
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the dimethylisoalloxazine core pattern; we use a simplified SMARTS string
    flavin_core_smarts = "Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C)c2cc1C"
    flavin_core = Chem.MolFromSmarts(flavin_core_smarts)
    
    # Check if the molecule contains the core structure
    if not mol.HasSubstructMatch(flavin_core):
        return False, "Does not contain the dimethylisoalloxazine core"
    
    # Get matches and check for substitution at the 10 position
    matches = mol.GetSubstructMatch(flavin_core)
    if matches:
        # Assuming position 10 corresponds to the nitrogen atom in the core
        position_10 = matches[5]  # example index for position 10, adjust as necessary
        atom = mol.GetAtomWithIdx(position_10)
        
        # Check for substitution at position 10 (it should have more than one bond)
        if len(atom.GetNeighbors()) <= 1:
            return False, "No substitution at position 10"
    
    return True, "Contains dimethylisoalloxazine core with substitution at position 10"

# Example usage
smiles = "Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C)c2cc1C"
result, reason = is_flavin(smiles)
print(result, reason)