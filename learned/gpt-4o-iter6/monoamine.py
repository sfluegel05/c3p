"""
Classifies: CHEBI:63534 monoamine
"""
from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine is defined as a compound which contains one amino group connected
    to an aromatic ring by a two-carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for the aromatic ring and amino group connected by two-carbon chain
    aromatic_ring_pattern = Chem.MolFromSmarts("c1ccccc1")
    two_carbon_amino_pattern = Chem.MolFromSmarts("[NH2]-C-C-c1ccccc1")

    # Check for an aromatic ring
    if not mol.HasSubstructMatch(aromatic_ring_pattern):
        return False, "No aromatic ring found"

    # Check for two-carbon chain with amino group connected to aromatic ring
    if not mol.HasSubstructMatch(two_carbon_amino_pattern):
        return False, "Missing two-carbon chain with amino group connected to aromatic ring"

    # If all checks pass, it is classified as a monoamine
    return True, "Contains amino group connected to aromatic ring by a two-carbon chain"