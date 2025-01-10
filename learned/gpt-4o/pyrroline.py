"""
Classifies: CHEBI:23763 pyrroline
"""
from rdkit import Chem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is defined as a dihydropyrrole, which is a five-membered ring
    containing a nitrogen atom with varying levels of saturation.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrroline, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES into RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns identifying pyrroline structures: five-membered rings with 1 nitrogen, allowing various substitutions
    pyrroline_patterns = [
        Chem.MolFromSmarts("C1CC=CN1"),  # Basic 2-pyrroline, simple unsaturated ring
        Chem.MolFromSmarts("C1C=CCN1"),  # Simple variant in the unsaturation position
        Chem.MolFromSmarts("C1CCC=CN1"), # Variants allowing full saturation, mimicking the dihydropyrrolic conditions
        Chem.MolFromSmarts("C1CN=CC1"),  # Allowing nitrogen-centered unsaturation
        Chem.MolFromSmarts("C1CCN=CC1"), # Including different nitrogen placement
        Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#7]~1"), # Generalized atomic pattern allowing carbon and one nitrogen in a ring
    ]

    # Ensure all SMARTS patterns are valid
    for idx, pattern in enumerate(pyrroline_patterns):
        if pattern is None:
            return False, f"Invalid SMARTS pattern at index {idx}"

    # Check against pyrroline structural motifs
    for pattern in pyrroline_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains the pyrroline structural motif (dihydropyrrole core)"
    
    return False, "Does not contain the pyrroline structural motif"

# Example testing with diverse compounds
test_smiles = [
    "C1CC=CN1",  # Core pyrroline structure
    "OC(=O)[C@@H]1CCC=N1",  # Slightly modified pyrroline structure
    "CCc1c(C)c2c([C@@H](C)C[C@@]2(C)O)c(C)c1C1=C\\C(C(=O)N1)=C1\\C=C(NC1=O)c1c(C)c2[C@H](C)C[C@](C)(O)c2c(C)c1CC"  # More complex structures
]

results = [is_pyrroline(smiles) for smiles in test_smiles]
for idx, result in enumerate(results):
    print(f"Test {idx + 1}: SMILES: {test_smiles[idx]} -> {result}")