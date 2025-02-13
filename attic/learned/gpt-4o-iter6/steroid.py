"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid contains the cyclopenta[a]phenanthrene carbon skeleton.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    # Cyclopenta[a]phenanthrene core pattern for steroids
    steroid_core_smarts = "C1CCC2C3CCC4=CC(=O)CC=C4C3CCC12"
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES format."

    # Check if the molecule contains the steroid core
    steroid_core_pattern = Chem.MolFromSmarts(steroid_core_smarts)
    
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Does not contain the cyclopenta[a]phenanthrene core structure."

    # Count atoms for methyl typically at C-10 and C-13
    n_methyl = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CH3]")))
    if n_methyl < 2:
        return False, "Insufficient methyl groups, typically found at C-10 and C-13."

    # Check for any additional typical alkyl group on C-17
    # This can be complex and requires detailed analysis
    alkyl_pattern = Chem.MolFromSmarts("[C](-[C&H2R4]-[C])")
    if mol.HasSubstructMatch(alkyl_pattern):
        return True, "Matches cyclopenta[a]phenanthrene core with alkyl groups typical of steroids."

    return False, "Matches core but lacks typical alkyl decorations or structure deviations typical in steroids."

# Some example SMILES for test
examples = [
    "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)[C@H](F)C[C@]34C)[C@@H]1CC[C@@H]2O",  # Example steroid
    "CC(=O)O[C@H]1CC[C@H]2[C@@H]3CC[C@H]4CCCC[C@]4(C)[C@H]3CC[C@]12C"         # Example steroid
]

for smiles in examples:
    result, reason = is_steroid(smiles)
    print(f"SMILES: {smiles} -> Is Steroid: {result}, Reason: {reason}")