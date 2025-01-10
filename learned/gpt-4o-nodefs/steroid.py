"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    Steroids are characterized by the cyclopentanoperhydrophenanthrene (CPPP)
    core structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the steroid core: cyclopentanoperhydrophenanthrene
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4CCCCC4C3C=C2C1")  # A typical steroid backbone

    # Check for the steroid pattern in the molecule
    if mol.HasSubstructMatch(steroid_pattern):
        return True, "Contains the cyclopentanoperhydrophenanthrene core structure"

    return False, "Does not match the steroid structural core"

# Example usage
example_smiles = "O1[C@@]23[C@]([C@@]4([C@]([C@]5([C@](CC3)([C@](Oc6ccccc6)(C5)C)[H])[H])(CC4)[H])[H])[H](C=O)CC3O)C(CC2C)C(=O)C=C1"
print(is_steroid(example_smiles))