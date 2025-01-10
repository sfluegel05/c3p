"""
Classifies: CHEBI:38958 indole alkaloid
"""
from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    It checks for the presence of an indole skeleton and other typical alkaloid features.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES into RDKit Molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for common indole pattern
    indole_patterns = [
        Chem.MolFromSmarts('c1c[nH]c2cccc2c1'),  # canonical indole
        Chem.MolFromSmarts('c1cnc2cccc2c1'),    # indole without hydrogen
        Chem.MolFromSmarts('c1c[n,nH]c2ccccc2c1')  # more flexible with options for N
    ]
    
    # Check if any indole patterns are present
    if not any(mol.HasSubstructMatch(pattern) for pattern in indole_patterns):
        return False, "No recognizable indole-like skeleton found"

    # Verify presence of additional nitrogen atoms; typical for complex alkaloids
    n_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if len(n_atoms) <= 1:
        return False, "Not enough nitrogen atoms typical of indole alkaloids"

    return True, "Molecule possesses an indole structure with characteristics of an alkaloid"

# Example use case for testing based on provided SMILES strings
example_smiles = "O=C1N2[C@H]([C@@]3(O[C@](C(N3[C@H]1CC(C)C)=O)(NC(=O)[C@@H]4C=C5C6=C7C(NC=C7C[C@H]5N(C4)C)=CC=C6)CC)O)CCC2"
result, reason = is_indole_alkaloid(example_smiles)
print(result, reason)