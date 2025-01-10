"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
from rdkit import Chem

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    A monoterpenoid indole alkaloid contains an indole ring and terpenoid-like features.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expanded SMARTS pattern for various indole rings
    # Include additional structures for indole-related features
    indole_patterns = [
        Chem.MolFromSmarts('c1c2[nH]c3c(c2ccc1)c[nH]c3'),  # Multi-fusion indole structures
        Chem.MolFromSmarts('c1cc2cc[nH]c2cc1'),  # Basic indole structure
        Chem.MolFromSmarts('c1ccc2nc3ccc4ccccc4c3[nH]c2c1'),  # Complex indole-like structures
        Chem.MolFromSmarts('c1c[nH]c2ccccc12'),  # Simple indole motif
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in indole_patterns):
        return False, "No indole ring found"

    # Refined SMARTS pattern for terpenoid features
    # Include common cyclic and acyclic patterns typical of monoterpenoids
    terpenoid_patterns = [
        Chem.MolFromSmarts('C1(C)C=CC=C1'),  # Acyclic monoterpenoid-like structure
        Chem.MolFromSmarts('C(C)(C)C=C'),  # Simple acyclic isoprene-like features
        Chem.MolFromSmarts('C1(C)CCC(C1)'),  # Monocyclic monoterpenoids
        Chem.MolFromSmarts('C1(CCCCC1)'),  # Saturated rings
        Chem.MolFromSmarts('C1CCC=C(C1)')  # Cyclic with double bonds
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in terpenoid_patterns):
        return False, "No terpenoid features found"

    # If molecule passes all checks, classify as a monoterpenoid indole alkaloid
    return True, "Contains both indole ring and terpenoid features."