"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
from rdkit import Chem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.
    An 11beta-hydroxy steroid contains a steroid core structure with a hydroxyl group at the 11beta position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Defining a general steroid core pattern
    steroid_core_pattern = Chem.MolFromSmarts("C1CCC2C3C(CCC4C3CCC2C1)C4")  # Generalized core

    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Steroid core structure not found"

    # Defining the 11-beta-hydroxy specific pattern, assuming stereochemistry is crucial here
    # Using a pattern to capture hydroxyl group on a steroid system
    hydroxy11b_pattern = Chem.MolFromSmarts("[C@H](O)[C@@H]1CCC2C(CCC3=C2[C@@H]1)C3")

    if not mol.HasSubstructMatch(hydroxy11b_pattern):
        return False, "Missing 11beta-hydroxy group"

    return True, "Contains steroid core structure with 11beta-hydroxy group"