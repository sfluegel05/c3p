"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    Steroids are characterized by the tetracyclic cyclopentanoperhydrophenanthrene (CPPP) core structure,
    which includes three six-membered rings and one five-membered ring.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a steroid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define various SMARTS patterns for the steroid core (CPPP structure)
    # This allows some flexibility in ring orientations
    steroid_patterns = [
        Chem.MolFromSmarts("C1CCC2C3CCC4CCCCC4C3C=C2C1"),  # Original pattern to be used as baseline
        Chem.MolFromSmarts("C1CC2C3CCC4C2CCC4C3C1"),  # Tetracyclic core without explicit double bonds
        Chem.MolFromSmarts("C1CC2CC3CCC4C(C)CCC4(C)C3C2C1"),  # Allowing substitutions and stereochemistry
        Chem.MolFromSmarts("C1CCC2CC(CC3CCC4CCCCC4C3)C2C1"),  # Variant pattern to account for different conformations
    ]

    # Check for any steroid pattern in the molecule
    for pattern in steroid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains the cyclopentanoperhydrophenanthrene core structure"

    return False, "Does not match the steroid structural core"

# Example usage
example_smiles = "O1[C@@]23[C@]([C@@]4([C@]([C@]5([C@](CC3)([C@](Oc6ccccc6)(C5)C)[H])[H])(CC4)[H])[H])[H](C=O)CC3O)C(CC2C)C(=O)C=C1"
print(is_steroid(example_smiles))