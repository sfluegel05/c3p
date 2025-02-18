"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: CHEBI:flavonols
Flavonols are hydroxyflavones with a hydroxyl group at position 3 of the heterocyclic (C) ring.
"""
from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol is a hydroxyflavone with a hydroxyl group at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for the flavone backbone with hydroxyl at position 3 of the C ring
    # The pattern matches a six-membered ring (pyrone) with:
    # - Oxygen at position 1
    # - Hydroxyl group at position 3
    # - Ketone at position 4
    flavonol_pattern = Chem.MolFromSmarts("[O;r1]1[C;r1][C;r1](O)[C;r1](=O)[C;r1][C;r1]1")

    # Check for the presence of the core structure
    if mol.HasSubstructMatch(flavonol_pattern):
        return True, "Contains flavone backbone with hydroxyl group at position 3 of the heterocyclic ring"
    else:
        return False, "Does not contain required flavone structure with hydroxyl at position 3"