"""
Classifies: CHEBI:23763 pyrroline
"""
from rdkit import Chem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is a dihydropyrrole with less unsaturation compared to pyrrole.

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

    # SMARTS pattern for a dihydropyrrole: five-membered ring with 1 nitrogen, with partial saturation
    pyrroline_pattern = Chem.MolFromSmarts("C1=CNCC1 | C1=CNC=C1 | C1=CNC(C)=C1")
    
    if mol.HasSubstructMatch(pyrroline_pattern):
        return True, "Contains the pyrroline structural motif (dihydropyrrole core)"
    else:
        return False, "Does not contain the pyrroline structural motif"

# Testing the function with given examples
test_smiles = [
    "C1CC=CN1",
    "OC(=O)[C@@H]1CCC=N1",
    "C1CC2=C(CC1)C(N(C2=O)C3=C(C=C(C(=C3)OC(C)C#C)Cl)F)=O",
    "CCc1c(C)c2c([C@@H](C)C[C@@]2(C)O)c(C)c1C1=C\\C(C(=O)N1)=C1\\C=C(NC1=O)c1c(C)c2[C@H](C)C[C@](C)(O)c2c(C)c1CC"
]

results = [is_pyrroline(smiles) for smiles in test_smiles]
for idx, result in enumerate(results):
    print(f"Test {idx + 1}: SMILES: {test_smiles[idx]} -> {result}")