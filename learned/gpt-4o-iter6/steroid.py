"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid generally consists of a tetracyclic structure with three six-membered
    rings and one five-membered ring (cyclopenta[a]phenanthrene skeleton).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a steroid, False otherwise
        str: Reason for classification
    """
    # General steroid backbone pattern with tolerance for natural variations
    steroid_core_smarts = "C1CC2(C)C3C(C)C4C=CC=C5=C4C3=C2C1=C5"  # Cyclopenta[a]phenanthrene core

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES format."

    # Check if the molecule contains the basic steroid core structure
    steroid_core_pattern = Chem.MolFromSmarts(steroid_core_smarts)
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Does not contain the cyclopenta[a]phenanthrene core structure."

    # Common functional groups found in steroids
    common_modifications = [
        Chem.MolFromSmarts("C=O"),  # Keto group
        Chem.MolFromSmarts("O"),    # Hydroxyl group or ether
        Chem.MolFromSmarts("C#O")   # Acetyl or other carbonyl-containing groups
    ]

    # Check for the presence of common functional modifications
    for feature in common_modifications:
        if mol.HasSubstructMatch(feature):
            return True, "Contains steroid core with characteristic functional groups."
    
    return False, "Contains the core structure but lacks any common functional groups indicative of steroids."

# Example usage with given examples
examples = [
    "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)[C@H](F)C[C@]34C)[C@@H]1CC[C@@H]2O",  # Steroid
    "CC(=O)O[C@H]1CC[C@H]2[C@@H]3CC[C@H]4CCCC[C@]4(C)[C@H]3CC[C@]12C"         # Steroid
]

for smiles in examples:
    result, reason = is_steroid(smiles)
    print(f"SMILES: {smiles} -> Is Steroid: {result}, Reason: {reason}")