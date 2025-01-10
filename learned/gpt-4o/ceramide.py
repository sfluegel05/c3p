"""
Classifies: CHEBI:17761 ceramide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.

    A ceramide consists of a sphingoid base (long unsaturated aliphatic chain with amino alcohol)
    and an amide-linked fatty acid, typically with a chain length from 14 to 26 carbons.

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

    # Look for sphingoid base-like structure: long-chain aliphatic amino alcohol
    sphingoid_pattern = Chem.MolFromSmarts("[C,C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][NH][C](O)")
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "No sphingoid base found"
    
    # Look for amide bond - linked fatty acid
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found"

    # Count carbon atoms in fatty acid chain - typically 14 to 26
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 14 or c_count > 26:
        return False, f"Fatty acid chain length {c_count} not in range (14-26)"

    return True, "Contains sphingoid base and amide-linked fatty acid"

# Test the function with given examples
examples = [
    "CCCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCCCC",  # N-hexacosanoylsphinganine
    "O1C(C(O)C(O)C(O)C1OCC(NC(=O)C(O)CCCCCCCCCCCC)C(O)/C=C\\CC/C=C\\CCCCCCCCC)CO",  # AS 1-1
    # Add more example SMILES strings to test
]

for smiles in examples:
    result, reason = is_ceramide(smiles)
    print(f"SMILES: {smiles} is a ceramide? {result}: {reason}")