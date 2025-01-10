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

    # Define the dimethylisoalloxazine core pattern with flexible bonding for the core structure
    dimethylisoalloxazine_core = Chem.MolFromSmarts("Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C)[c]2cc1C")

    # Check if the molecule contains the core structure
    core_match = mol.GetSubstructMatch(dimethylisoalloxazine_core)
    if not core_match:
        return False, "Does not contain the dimethylisoalloxazine core"

    # Ensure there is a substitution at the position bonded off the recognized core structure
    # Find potential attachment points from the substitution view
    core_atoms = set(core_match)
    for atom in mol.GetAtoms():
        if atom.GetIdx() in core_atoms and atom.GetSymbol() == 'C':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in core_atoms:
                    # Check if it's directly connected to the core but not part of the core
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