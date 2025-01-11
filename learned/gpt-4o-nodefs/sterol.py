"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol typically has a fused four-ring structure (cyclopentanoperhydrophenanthrene) with a hydroxyl group on the C3 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core sterol structure using a SMARTS pattern
    sterol_pattern = Chem.MolFromSmarts("C1(C(C2CCC3C(C(C(C4=CC=C(C(C34)C2)C1)C)O)O)C)CCCC(C)C")

    # Check if the molecule contains the sterol pattern
    if mol.HasSubstructMatch(sterol_pattern):
        return True, "Contains the characteristic four-ring structure of a sterol"

    return False, "Does not contain the characteristic sterol structure"