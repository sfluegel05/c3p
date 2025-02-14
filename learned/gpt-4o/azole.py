"""
Classifies: CHEBI:68452 azole
"""
from rdkit import Chem

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole based on its SMILES string.
    An azole is defined as a monocyclic heteroarene with a five-membered ring
    containing at least one nitrogen atom and possibly other heteroatoms such
    as nitrogen, sulfur, or oxygen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an azole, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for different azoles
    azole_patterns = [
        Chem.MolFromSmarts("c1ncnc1"),  # Pyrazole
        Chem.MolFromSmarts("c1ncc[nH]1"),  # Imidazole
        Chem.MolFromSmarts("c1nncn1"),  # 1,2,3-Triazole
        Chem.MolFromSmarts("c1nnnc1"),  # Tetrazole
        Chem.MolFromSmarts("c1ncn[s,o]1"),  # Thiazole/Oxazole
    ]

    # Check if any azole pattern is present in the molecule
    for pattern in azole_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains an azole ring based on defined patterns"

    return False, "No azole ring found"