"""
Classifies: CHEBI:16385 organic sulfide
"""
from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide based on its SMILES string.
    An organic sulfide has a sulfur atom bonded to two organic groups, designated as R1-S-R2,
    where neither R1 nor R2 is a hydrogen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic sulfide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern: Match sulfur with likely organic surroundings
    # This pattern considers more versatile scenarios including complex organic attachments
    sulfide_patterns = [
        "[C,c][S][C,c]",  # Basic organic paths, allows aromatic or non-aromatic
        "[a][S][!#1]",    # Aromatic atoms connected to sulfur, not hydrogen
        "[!#1][S][!#1]"   # More general check where sulfur is not bonded to any hydrogen
    ]

    # Check molecule for organic sulfide structure using multiple patterns
    for pattern in sulfide_patterns:
        sulfide_match = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(sulfide_match):
            return True, f"Contains organic sulfide group matching pattern: {pattern}"

    return False, "No organic sulfide group found"

# Test with various SMILES examples
test_smiles = [
    "CC(=O)N[C@@H]1[C@@H](O)[C@H](O[C@@H]2O[C@H](CO)[C@H](O[C@H]3O[C@H](CO)[C@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@@H](CO)O[C@H]1OCCSCCC(N)=O",
    "C\\C=C\\S\\C=C\\C",
    "CSC1=N[C@](C)(C(=O)N1Nc1ccccc1)c1ccccc1",
    "S(C(C)C)C(C)C",
    "CC1=C(C=C(C=C1)SC2=C(C=CC=N2)C(=O)O)C",
    "CCOP(=S)(OCC)SCSc1ccc(Cl)cc1",
]

for smiles in test_smiles:
    result, reason = is_organic_sulfide(smiles)
    print(f"SMILES: {smiles} -> Is Organic Sulfide: {result}, Reason: {reason}")