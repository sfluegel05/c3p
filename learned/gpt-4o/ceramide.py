"""
Classifies: CHEBI:17761 ceramide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    Ceramides are characterized by a sphingoid base with an amide-linked fatty acid.
    The sphingoid base typically has an amino alcohol, and the fatty acids tend to range
    from 14 to 26 carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ceramide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sphingoid base structure - Long aliphatic chain with aminodiol pattern
    sphingoid_pattern = Chem.MolFromSmarts("C[C@@H](O)[C@@H](O)CO")
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "No sphingoid base found"

    # Look for an amide linkage pattern
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found"

    # Calculate the total number of carbon atoms
    carbon_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check if total carbon count meets expectations for ceramides
    if carbon_atoms < 18:
        return False, f"Too few carbon atoms for ceramide, found {carbon_atoms}"

    return True, "Contains sphingoid base with amide-linked fatty acid"

# Test the function
examples = [
    "CCCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCCCC",  # N-hexacosanoylsphinganine
    "O1C(C(O)C(O)C(O)C1OCC(NC(=O)C(O)CCCCCCCCCCCC)C(O)/C=C\\CC/C=C\\CCCCCCCCC)CO",  # AS 1-1
    # Add more example SMILES strings to test
]

for smiles in examples:
    result, reason = is_ceramide(smiles)
    print(f"SMILES: {smiles} is a ceramide? {result}: {reason}")