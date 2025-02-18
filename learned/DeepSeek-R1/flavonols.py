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

    # Corrected SMARTS pattern for flavonol core structure:
    # Six-membered heterocyclic ring (C ring) with:
    # - Oxygen atom in the ring (position 1)
    # - Hydroxyl group at position 3
    # - Ketone at position 4
    # All atoms part of a 6-membered ring (r6)
    flavonol_pattern = Chem.MolFromSmarts("[O;r6]1-[C;r6]-[C;r6](O)-[C;r6](=O)-[C;r6]-[C;r6]-1")

    # Check for the presence of the core structure
    if mol.HasSubstructMatch(flavonol_pattern):
        return True, "Contains flavone backbone with hydroxyl group at position 3 of the heterocyclic ring"
    else:
        return False, "Does not contain required flavone structure with hydroxyl at position 3"