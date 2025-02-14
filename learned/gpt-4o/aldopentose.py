"""
Classifies: CHEBI:33916 aldopentose
"""
from rdkit import Chem

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    An aldopentose is a pentose with a terminal aldehyde group or existing in
    a cyclic (furanose or pyranose) form characteristic of aldopentoses, excluding
    any substitutions that alter the core sugar structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldopentose, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms, excluding any attached to phosphate (P) or other major substituent groups
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:
        return False, f"Expected at least 5 carbons, found {c_count}"

    # Exclude structures with phosphate or major substituent groups
    unwanted_patterns = [
        Chem.MolFromSmarts("OP(=O)(O)"),  # Phosphate groups
        Chem.MolFromSmarts("C(=O)[O-]"),  # Carboxylate, found in uronic acids
    ]
    if any(mol.HasSubstructMatch(pattern) for pattern in unwanted_patterns):
        return False, "Contains unwanted substituent groups like phosphate"

    # Check for linear form with terminal aldehyde group (-CHO)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[C]")
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Contains a terminal aldehyde group applied to a pentose core"

    # Detect potential pyranose form (six-membered ring with 5 carbons and 1 oxygen)
    pyranose_pattern = Chem.MolFromSmarts("C1OC[C@@H](O)[C@H](O)[C@H]1O")  # Better encapsulate the hydroxyl pattern
    if mol.HasSubstructMatch(pyranose_pattern):
        return True, "Contains a pyranose ring structure"

    # Detect potential furanose form (five-membered ring with 4 carbons, 1 oxygen)
    furanose_pattern = Chem.MolFromSmarts("C1C[C@@H](O[C@@H]1O)O")
    if mol.HasSubstructMatch(furanose_pattern):
        return True, "Contains a furanose ring structure"

    return False, "Does not match aldopentose profile"