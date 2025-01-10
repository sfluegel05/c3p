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

    # More detailed steroid core pattern
    steroid_core_patterns = [
        "C[C@H]1CC[C@H]2[C@@H]3CC[C@H]4C(=O)C=C[C@@]4(C)[C@@H]3CC[C@]12C",  # Fully stereospecific pattern
        "C1CCC2C1CCC3C2CCC4C3CCC4",  # Typical steroid stereochemistry
        "C[C@@H]1CC[C@@H]2[C@@H]3CC[C@H]4C(=O)OCCCC4(C)[C@@H]3CCC2C1"  # Alternate forms
    ]
    if not any(mol.HasSubstructMatch(Chem.MolFromSmarts(p)) for p in steroid_core_patterns):
        return False, "Steroid core structure not found"

    # 11beta-hydroxy group: targeting specific position in the steroid skeletal
    hydroxyl_pattern = Chem.MolFromSmarts("[C@]([C@H])([OH])")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing 11beta-hydroxy group"

    return True, "Contains steroid core structure with 11beta-hydroxy group"