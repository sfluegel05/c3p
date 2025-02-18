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

    # Define the dimethylisoalloxazine core pattern; Start by capturing variations
    flavin_core_smarts = "Cc1cc2nc3c(nc(=O)[nH]c3=O)n(c2cc1[C,N,O,S])[C,N,O,S]"  # Adjusted for more flexible core matching
    flavin_core = Chem.MolFromSmarts(flavin_core_smarts)
    
    # Check if the molecule contains the core structure
    if not mol.HasSubstructMatch(flavin_core):
        return False, "Does not contain the dimethylisoalloxazine core"
    
    # Get matches for the core structure
    matches = mol.GetSubstructMatch(flavin_core)
    if matches:
        # Check for substitution at the position typically bonded to the core at position 10
        core_atoms = set(matches)
        for atom in mol.GetAtoms():
            if atom.GetIdx() not in core_atoms and atom.GetSymbol() in ['C', 'N', 'O', 'S']:  # Ensure substitution involves valid atoms
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() in core_atoms:
                        return True, "Contains dimethylisoalloxazine core with substitution at position 10"
    
    return False, "No substitution at position 10 of the dimethylisoalloxazine core"

# Example usage
examples = [
    "Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C)c2cc1C",  # Known base flavin structure
    "Cc1cc2N(C3=NC(=O)NC(=O)C3=N2)C[C@H](O)[C@H](O)[C@H](O)COC4OC(CO)C(C4O)O"  # Example of flavin with complex substitution
]
for smiles in examples:
    result, reason = is_flavin(smiles)
    print(f"SMILES: {smiles}\nResult: {result}, Reason: {reason}\n")