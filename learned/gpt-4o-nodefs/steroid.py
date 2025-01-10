"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    Steroids are characterized by the tetracyclic cyclopentanoperhydrophenanthrene (CPPP) core structure,
    which includes three six-membered rings and one five-membered ring in various conformations.
    
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
    # This includes several common variations to allow for different steroid conformations
    steroid_patterns = [
        Chem.MolFromSmarts("C1CCC2C3CCC4C(C)CCC4C3C=C2C1"),  # Basic CPPP pattern
        Chem.MolFromSmarts("C1CCC2C3CCC4C(C)CCC4C3CC2C1"),  # Variant with different orientations
        Chem.MolFromSmarts("C1CC2CC3CCC4C(C)CCC4(C)C3C2C1"),  # Allowing some random substitutions
        Chem.MolFromSmarts("C1CC2C3(C)CCC4C(CC)(C)CCC4(CC3C=C2)C1"),  # Steroid-like flexibility
    ]
    
    # Check for any steroid pattern in the molecule
    for pattern in steroid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains the cyclopentanoperhydrophenanthrene core structure indicative of steroids"
    
    return False, "Does not match the steroid structural core"

# Example usage
example_smiles = "O=C1[C@@]2(C(=CC[C@H]1[C@H]3[C@@]([C@@H]([C@@H](CC(=O)C(C(C)C)C)C)CC3)(CCO)C)C[C@@H](O)CC2)C"
print(is_steroid(example_smiles))