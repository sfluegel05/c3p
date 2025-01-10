"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    Steroids have a cyclopenta[a]phenanthrene carbon skeleton, with various modifications.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a steroid, False otherwise
        str: Reason for classification
    """
    # SMARTS pattern for the core cyclopenta[a]phenanthrene skeleton
    steroid_core_smarts = "[#6]1[#6][#6]2[#6]3([#6](C)C4CC5C3C2C[C@]14)CC5" 

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES format."

    # Check if the molecule matches the steroid core pattern
    steroid_core_pattern = Chem.MolFromSmarts(steroid_core_smarts)
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Does not contain the cyclopenta[a]phenanthrene core structure."

    # Optional: Presence of functional groups that are characteristic of steroids
    common_modifications = [
        Chem.MolFromSmarts("[C@H](O)"),
        Chem.MolFromSmarts("C(=O)"),
        Chem.MolFromSmarts("C[C@H](O)C"),
        Chem.MolFromSmarts("CC(C)"),
        Chem.MolFromSmarts("C(=O)O"),
    ]

    # Check for at least one common steroid modification
    for feature in common_modifications:
        if mol.HasSubstructMatch(feature):
            return True, "Matches the core with common steroid modifications."
    
    return True, "Matches the core steroid structure, but lacks common functional modifications."

# Example usage with some example SMILES strings
examples = [
    "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)[C@H](F)C[C@]34C)[C@@H]1CC[C@@H]2O",  # Steroid
    "CC(=O)O[C@H]1CC[C@H]2[C@@H]3CC[C@H]4CCCC[C@]4(C)[C@H]3CC[C@]12C"         # Steroid
]

for smiles in examples:
    result, reason = is_steroid(smiles)
    print(f"SMILES: {smiles} -> Is Steroid: {result}, Reason: {reason}")