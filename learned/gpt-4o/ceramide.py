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

    # Revised pattern: Long-chain base (aliphatic chain, aminodiol)
    sphingoid_pattern = Chem.MolFromSmarts("C[C@@H](O)[C@@H](O)CO[C@H]([NH,C](C(=O)))")
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "No sphingoid base found"

    # Look for amide linkage - linked fatty acid
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found"

    # Assess chain length - count carbon atoms in associated groups
    carbon_sn1 = rdMolDescriptors.CalcNumAliphaticCarbocycles(mol)
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if total_carbons < 18:  # Minimal total carbons to account for backbone plus fatty acid chain
        return False, f"Too few carbon atoms, found {total_carbons}"

    return True, "Contains sphingoid base with amide-linked fatty acid"

# Test the function
examples = [
    "CCCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCCCC",  # N-hexacosanoylsphinganine
    # Add more example SMILES strings to test
]

for smiles in examples:
    result, reason = is_ceramide(smiles)
    print(f"SMILES: {smiles} is a ceramide? {result}: {reason}")