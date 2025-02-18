"""
Classifies: CHEBI:16158 steroid sulfate
"""
"""
Classifies: CHEBI:37962 steroid sulfate
"""
from rdkit import Chem

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    A steroid sulfate is a sulfuric ester formed by the condensation of a steroid's hydroxyl group with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid sulfate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sulfate ester group (-OSO3)
    sulfate_pattern = Chem.MolFromSmarts("[OX2]S(=O)(=O)[OX2]")
    if not mol.HasSubstructMatch(sulfate_pattern):
        return False, "No sulfate ester group detected"

    # Check for steroid nucleus (four fused rings)
    # SMARTS pattern for the steroid skeleton (simplified)
    steroid_pattern = Chem.MolFromSmarts(
        "[C]1CC[C@]2[C@]1([C@H]([C@H]3[C@H]([C@@H]4[C@@](CC3)([C@H](CC4)CC2)C)C)C)C"
    )
    if not mol.HasSubstructMatch(steroid_pattern):
        # Try alternative pattern without stereochemistry
        steroid_pattern_general = Chem.MolFromSmarts(
            "C12CCC3C4CC(CC(C4)C3)C1CCC2"
        )
        if not mol.HasSubstructMatch(steroid_pattern_general):
            return False, "No steroid nucleus detected"

    return True, "Contains steroid nucleus with sulfate ester group"