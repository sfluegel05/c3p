"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid contains the cyclopenta[a]phenanthrene carbon skeleton with possible
    extensions for methyl groups at C-10 and C-13, and often alkyl groups at C-17.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    # Steroids have a cyclopenta[a]phenanthrene core
    steroid_core_smarts = "C1CC[C@H]2C[C@H]3C[C@H](CC4C3CCC4)[C@H]2C1"  # General core with variations

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES format."

    # Check if the molecule contains the core structure that can undergo modifications
    steroid_core_pattern = Chem.MolFromSmarts(steroid_core_smarts)
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Does not contain the cyclopenta[a]phenanthrene core structure."

    # Check for methyl groups usually at C-10 and C-13
    methyl_pattern = Chem.MolFromSmarts("[C]([C&H3])[C&H]")  # Simplified pattern for methyl groups on ring system
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)
    if len(methyl_matches) < 2:
        return False, "Insufficient methyl groups, typically found at C-10 and C-13."

    # Alkyl group pattern often at C-17
    alkyl_group_pattern = Chem.MolFromSmarts("[C]([C&H2R4]-[C])")  # Alkyl group connected to ring
    if mol.HasSubstructMatch(alkyl_group_pattern):
        return True, "Matches cyclopenta[a]phenanthrene core with alkyl groups typical of steroids."

    return False, "Matches core but lacks typical alkyl decorations or substitutions seen in steroids."

# Some example SMILES for test
examples = [
    "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)[C@H](F)C[C@]34C)[C@@H]1CC[C@@H]2O",  # Steroid
    "CC(=O)O[C@H]1CC[C@H]2[C@@H]3CC[C@H]4CCCC[C@]4(C)[C@H]3CC[C@]12C"         # Steroid
]

for smiles in examples:
    result, reason = is_steroid(smiles)
    print(f"SMILES: {smiles} -> Is Steroid: {result}, Reason: {reason}")