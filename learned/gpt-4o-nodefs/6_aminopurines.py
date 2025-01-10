"""
Classifies: CHEBI:20706 6-aminopurines
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule belongs to the 6-aminopurines class based on its SMILES string.
    6-aminopurines are characterized by a purine ring with an amino group at the 6-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) - True if molecule is a 6-aminopurine, False otherwise with reason
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the 6-aminopurine SMARTS pattern
    # A purine with an amino group at the 6-position, ignoring potential tautomers and substituents
    aminopurine_pattern = Chem.MolFromSmarts("Nc1ncnc2c1ncnc2")

    # Check for 6-aminopurine substructure
    if mol.HasSubstructMatch(aminopurine_pattern):
        return True, "Contains 6-aminopurine substructure"

    return False, "Does not contain 6-aminopurine substructure"