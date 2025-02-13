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
        Chem.MolFromSmarts("C1CC=CN1"),  # 2,3-dihydro-1H-pyrrole
        Chem.MolFromSmarts("C1C=CCN1"),  # 1,2-dihydro-1H-pyrrole
        Chem.MolFromSmarts("C1CN=CC1"),  # 2,5-dihydro-1H-pyrrole
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

# Testing the function with given examples
test_smiles = [
    "C1CC=CN1",  # Positive case
    "OC(=O)[C@@H]1CCC=N1",  # Positive case
    "C1CC2=C(CC1)C(N(C2=O)C3=C(C=C(C(=C3)OC(C)C#C)Cl)F)=O",  # Potentially positive case
    "CCc1c(C)c2c([C@@H](C)C[C@@]2(C)O)c(C)c1C1=C\\C(C(=O)N1)=C1\\C=C(NC1=O)c1c(C)c2[C@H](C)C[C@](C)(O)c2c(C)c1CC"  # Potentially negative or complex
]

results = [is_pyrroline(smiles) for smiles in test_smiles]
for idx, result in enumerate(results):
    print(f"Test {idx + 1}: SMILES: {test_smiles[idx]} -> {result}")