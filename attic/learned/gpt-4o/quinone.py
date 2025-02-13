"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    Quinones are characterized by having a fully conjugated cyclic dione structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define quinone patterns
    quinone_patterns = [
        Chem.MolFromSmarts("O=C1C=CC(=O)C=C1"),  # Benzoquinone core
        Chem.MolFromSmarts("O=C1C=CC(=O)c2ccccc12"),  # Naphthoquinone core
        Chem.MolFromSmarts("O=C1C=CC(=O)c2c1cccc2"),  # Anthraquinone core
        Chem.MolFromSmarts("O=C1c2ccccc2C(=O)c3ccccc13"),  # Phenanthrenequinone
    ]

    # Check for matches
    for pattern in quinone_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Found quinone pattern"

    return False, "No quinone pattern found"