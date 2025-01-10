"""
Classifies: CHEBI:23763 pyrroline
"""
from rdkit import Chem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is a dihydropyrrole with less unsaturation compared to pyrrole,
    consisting of a five-membered ring with a nitrogen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrroline, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns for pyrroline: five-membered ring with 1 nitrogen and varying saturation
    pyrroline_patterns = [
        Chem.MolFromSmarts("[NH]1CCC=C1"),  # 1H-pyrrolines, considering basic saturation
        Chem.MolFromSmarts("[NH]1C=CCC1"),  # Another possible hydrogenation variant
        Chem.MolFromSmarts("[NH]1CC=CC1"),  # 2,5-dihydro-1H-pyrrole structure
        Chem.MolFromSmarts("[nH]1CCC=C1"),  # Higher unsaturation with aromatic N
        Chem.MolFromSmarts("[nH]1C=CC=C1"), # Additional variant
    ]

    # Validate each SMARTS compilation
    for idx, pattern in enumerate(pyrroline_patterns):
        if pattern is None:
            return False, f"Invalid SMARTS pattern at index {idx}"

    # Check for any matching pyrroline structural motif
    for pattern in pyrroline_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains the pyrroline structural motif (dihydropyrrole core)"
    
    return False, "Does not contain the pyrroline structural motif"

# Example testing with edge cases
test_smiles = [
    "C1CC=CN1",  # Core pyrroline structure
    "OC(=O)[C@@H]1CCC=N1",  # Modified pyrroline structure
    "CCc1c(C)c2c([C@@H](C)C[C@@]2(C)O)c(C)c1C1=C\\C(C(=O)N1)=C1\\C=C(NC1=O)c1c(C)c2[C@H](C)C[C@](C)(O)c2c(C)c1CC"  # Complicated structure
]

results = [is_pyrroline(smiles) for smiles in test_smiles]
for idx, result in enumerate(results):
    print(f"Test {idx + 1}: SMILES: {test_smiles[idx]} -> {result}")