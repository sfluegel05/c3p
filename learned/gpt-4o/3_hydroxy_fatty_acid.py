"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid has a hydroxyl group at the beta- or 3-position relative to
    the carboxyl group and is part of a long hydrocarbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for 3-hydroxy fatty acids
    # Carboxylic acid group at one end: C(=O)O
    # Hydroxyl group at the 3rd position from the carboxyl group: [C@H](O) or [C@@H](O) indicates potential stereochemistry
    pattern = Chem.MolFromSmarts("C(=O)O-[CH2]-[CH](O)-[CH2]-[CH3]")

    # Check for the pattern match
    if mol.HasSubstructMatch(pattern):
        return True, "Contains hydroxyl group at the 3-position and carboxylic acid group"

    # Additional check for carbon chain length > 10 to ensure fatty acid characteristic
    carbon_chain_length = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_chain_length <= 10:
        return False, f"Carbon chain too short for fatty acid, got length: {carbon_chain_length}"

    return False, "Does not exhibit necessary hydroxyl and carboxyl pattern at the correct position"