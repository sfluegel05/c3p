"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    Sesquiterpenoids are characterized by a C15 backbone that may be rearranged or modified.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Adjust the expected carbon count range to accommodate structural rearrangements.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 13 or c_count > 20:  # Wider range to include derivatives
        return False, f"Expected around 15 carbons, got {c_count}"

    # Refined set of structural and functional patterns that may indicate sesquiterpenoids
    substructures = [
        Chem.MolFromSmarts('C1CCC(C)C1'),  # Minimal cycloalkane ring
        Chem.MolFromSmarts('O=C'),  # Carbonyl groups
        Chem.MolFromSmarts('C=C'),  # Double bonds
        Chem.MolFromSmarts('[OH]'),  # Hydroxy group
        Chem.MolFromSmarts('C1=CC=CC=C1'),  # Aromatic ring possible in complex sesquiterpenoids
        Chem.MolFromSmarts('O=C(O)'),  # Carboxylic acid or ester moiety
        Chem.MolFromSmarts('C1=C[C@H](C)C[C@@H]1'),  # Complex stereo substructure
    ]

    # Check for several sesquiterpenoid characteristrics
    matches = sum(mol.HasSubstructMatch(pattern) for pattern in substructures)

    # Empirically decide how many matches define a sesquiterpenoid
    if matches >= 3:
        return True, "Contains multiple sesquiterpenoid characteristic features"

    return False, "Does not exhibit enough sesquiterpenoid characteristic features"