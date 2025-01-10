"""
Classifies: CHEBI:17761 ceramide
"""
from rdkit import Chem

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    Ceramides are characterized by a sphingoid base with an amide-linked fatty acid.
    The sphingoid base typically has an amino alcohol, and the fatty acids tend to range
    from 14 to 26 carbons. The structure often includes a hydroxyl group on carbon 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ceramide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Look for a long aliphatic chain with a terminal NH group, suggesting a sphingoid base
    sphingoid_pattern = Chem.MolFromSmarts("N[C@H](CO)[C@H](O)CC")  # Generic for long chain amino alcohol
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "No sphingoid base found"

    # Look for an amide linkage pattern
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found"

    # Calculate the total number of carbon atoms
    carbon_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Count carbon atoms in the fatty acid chain portion to ensure it's within ceramide range
    fatty_acid_carbon_count = sum(1 for atom in mol.GetAtoms()
                                  if atom.GetAtomicNum() == 6 and atom.GetIsAromatic() == False)

    # Check if total fatty acid carbon count is within the expected range for ceramides
    if not (14 <= fatty_acid_carbon_count <= 26):
        return False, f"Fatty acid chain has {fatty_acid_carbon_count} carbons, outside ceramide range"

    # Check for presence of a commonly observed hydroxyl group
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if hydroxyl_count < 2:
        return False, "Hydroxyl groups found less than expected for ceramide"

    return True, "Contains sphingoid base with amide-linked fatty acid"

# Test the function
examples = [
    "CCCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCCCC",  # N-hexacosanoylsphinganine
    "O1C(C(O)C(O)C(O)C1OCC(NC(=O)C(O)CCCCCCCCCCCC)C(O)/C=C\\CC/C=C\\CCCCCCCCC)CO",  # AS 1-1
    "CCCCCCCCCCCCCCCCCC(=O)C1=CC=C(C=C1)NC(CCCCCCCCCCCCCCC)C(O)CC",
    "CCCCCCCCCCCCCCCCCCCCCCCCCCC=O"  # Control test not ceramide
]

for smiles in examples:
    result, reason = is_ceramide(smiles)
    print(f"SMILES: {smiles} is a ceramide? {result}: {reason}")