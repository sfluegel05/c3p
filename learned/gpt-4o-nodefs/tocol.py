"""
Classifies: CHEBI:39437 tocol
"""
from rdkit import Chem

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.
    Tocols such as tocopherols and tocotrienols have a chroman/chromene structure, long hydrocarbon tail, and hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tocol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Adjusted pattern for chroman/chromene core
    chroman_core_pattern = Chem.MolFromSmarts("C1[C@@H]2CCC[C@@]10CC=C(C(OC0=C2C)C)C")
    if not mol.HasSubstructMatch(chroman_core_pattern):
        return False, "No discernible chroman or chromene moiety found"

    # Check for long isoprenoid side chains with repetitive methyl substituents, common in tocotrienols
    long_hydrocarbon_pattern = Chem.MolFromSmarts("C(C)C(C)C")  # Match repeated methyl groups typical of isoprenoid tails
    tail_matches = mol.GetSubstructMatches(long_hydrocarbon_pattern)
    if len(tail_matches) < 2:
        return False, "Couldn't identify characteristic long isoprenoid tail with repetitive methyls"

    # Ensure at least one hydroxyl group attached to the aromatic core
    hydroxyl_pattern = Chem.MolFromSmarts("c(O)c")  # Hydroxyl group on aromatic carbon
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing aromatic hydroxyl group connection"

    return True, "Contains chroman/chromene moiety with key hydrocarbon tail and hydroxyl group"

# Testing for some cases
example_smiles = [
    "CC(C)CCC[C@H](C)CCC[C@H](C)CCC[C@]1(C)CCC2=C(C)C(O)=C(C)C(C)=C2O1",  # Example of valid tocol
    "ClC1=C2C=CC3=C1C4=CC=C[C@]4(O[C@@H]5OC([C@@H](N(C)C)[C@H]([C@H]5O)O)(C)C)[C@@H]3OC6=C(Cl)C=C([C@H](CC(OC[C@@H]2O)=O)N)C=C6O"  # Example of invalid tocol
]

# Check examples
for smiles in example_smiles:
    result, reason = is_tocol(smiles)
    print(f"SMILES: {smiles} -> Is Tocol: {result}, Reason: {reason}")