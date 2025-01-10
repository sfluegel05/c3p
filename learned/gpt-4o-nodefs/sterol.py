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

    # Define a more generalized core sterol structure using a SMARTS pattern
    # - Four-ring steroid backbones are captured
    # - Hydroxyl group at position 3
    sterol_pattern = Chem.MolFromSmarts("[#6]1[#6][#6]2[#6]3[#6]([#6]14)[#6][#6][C@H](O)[#6]3[#6]([#6]2)[#6]4")

    # Check if the molecule contains the sterol pattern
    if mol.HasSubstructMatch(sterol_pattern):
        return True, "Contains the characteristic sterol structure with four rings and a hydroxyl group"

    return False, "Does not contain the characteristic sterol structure"