"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid contains a modified cyclopenta[a]phenanthrene carbon skeleton.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    # Cyclopenta[a]phenanthrene core SMARTS with flexibility for common variations
    steroid_core_smarts = "C1CCC2C3CCC4CCCC[C@H]4[C@H]3[C@H]2C1"  # Classic representation of the backbone
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES format."

    # Check if the molecule contains the flexible steroid core structure
    steroid_core_pattern = Chem.MolFromSmarts(steroid_core_smarts)
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Does not contain the cyclopenta[a]phenanthrene core structure."

    # Further checks for variations in the core structure, such as methyl and alkyl groups
    # Methyl groups and other features can exist in diverse configurations
    extended_structure_patterns = [
        Chem.MolFromSmarts("CC[C@]1(C)C[C@H](CC2=C1CCC3C2CCC4[C@@H]3CC[C@@H]4C)")
    ]
    # Check for each extended structure
    for pattern in extended_structure_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains characteristic substructures of steroids including extensions."

    return True, "Matches cyclopenta[a]phenanthrene core."

# Some example SMILES for test
examples = [
    "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)[C@H](F)C[C@]34C)[C@@H]1CC[C@@H]2O",  # Steroid
    "CC(=O)O[C@H]1CC[C@H]2[C@@H]3CC[C@H]4CCCC[C@]4(C)[C@H]3CC[C@]12C"         # Steroid
]

for smiles in examples:
    result, reason = is_steroid(smiles)
    print(f"SMILES: {smiles} -> Is Steroid: {result}, Reason: {reason}")