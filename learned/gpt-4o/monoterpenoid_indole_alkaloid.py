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

    # Broader SMARTS pattern for various indole ring substructures
    indole_patterns = [
        Chem.MolFromSmarts('c1cc2c[nH]c2cc1'),  # Basic indole structure
        Chem.MolFromSmarts('c1ccc2c(c1)[nH]c3c2cccc3')  # Expanded indole-like structures
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in indole_patterns):
        return False, "No indole ring found"

    # SMARTS pattern for terpenoid-like structures (complex isoprene-related units)
    terpenoid_patterns = [
        Chem.MolFromSmarts('C1(C)C=CC=C1'),  # Acyclic monoterpenoid-like structure
        Chem.MolFromSmarts('C(C)(C)C=C'),  # More generic pattern already tried
        Chem.MolFromSmarts('C1(C)C2=CCC=C2C1')  # Cyclic terpenoid rings
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in terpenoid_patterns):
        return False, "No terpenoid features found (missing complex isoprene-related structures)"

    # If molecule passes all checks, classify as a monoterpenoid indole alkaloid
    return True, "Contains both indole ring and terpenoid features."