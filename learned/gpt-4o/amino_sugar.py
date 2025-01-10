"""
Classifies: CHEBI:28963 amino sugar
"""
from rdkit import Chem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is characterized by one or more amino groups replacing hydroxyl groups on a sugar molecule.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for sugar-like structure and amino group
    # A simple hexose-like sugar backbone pattern (C-C-C...O-O) can be used as a base
    # Note: This SMARTS pattern is generic and may not exactly match all sugar-like structures
    sugar_pattern = Chem.MolFromSmarts("C1OC[C@@H]([C@H](O)C1)O")  # llustrative sugar backbones, can be refined
    amino_group_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")  # Primary or secondary amine not in amides

    # Check for sugar-like backbone
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar backbone found"

    # Find all amino groups
    amino_matches = mol.GetSubstructMatches(amino_group_pattern)

    if not amino_matches:
        return False, "No amino groups replacing hydroxyl groups"

    # Further checks for explicit hydroxy displacement may include analyzing atom bonds directly
    # This requires investigating the specific bonding environment in the molecular structure

    return True, "Contains sugar backbone with one or more amino groups replacing hydroxy groups"