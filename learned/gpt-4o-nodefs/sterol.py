"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol typically has a fused four-ring structure (cyclopentanoperhydrophenanthrene),
    usually with a hydroxyl group at position 3 and a variable aliphatic side chain.

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

    # Define a core sterol structure using a more inclusive SMARTS pattern
    sterol_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4(C)C(C3C2C1)CC[C@H](O)C4")  # more representative sterol pattern

    # Check if the molecule contains the sterol pattern
    if mol.HasSubstructMatch(sterol_pattern):
        return True, "Contains the characteristic sterol structure with four rings and a hydroxyl group"

    return False, "Does not contain the characteristic sterol structure"