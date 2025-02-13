"""
Classifies: CHEBI:35358 sulfonamide
"""
from rdkit import Chem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide is an amide of a sulfonic acid RS(=O)2NR'2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfonamide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define sulfonamide pattern: RS(=O)2NR'
    sulfonamide_pattern = Chem.MolFromSmarts("S(=O)(=O)N")

    # Check for the sulfonamide substructure
    if mol.HasSubstructMatch(sulfonamide_pattern):
        return True, "Contains sulfonamide substructure"
    else:
        return False, "Does not contain sulfonamide substructure"

# Example usage
example_smiles = "C1CCC(CC1)NC(=O)C[C@@H]2C[C@H]3[C@@H]([C@@H](O2)CO)OC4=C3C=C(C=C4)NS(=O)(=O)C5=CC=C(C=C5)F"
result, reason = is_sulfonamide(example_smiles)
print(f"Is sulfonamide: {result}, Reason: {reason}")