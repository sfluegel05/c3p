"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    Quinones are characterized by a fully conjugated cyclic dione structure,
    derived from aromatic compounds by converting an even number of -CH= groups
    into -C(=O)- groups.

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
    
    # SMARTS patterns for various quinone structures
    quinone_patterns = [
        Chem.MolFromSmarts("c1c(=O)c(=O)c1"),  # Simple benzoquinone, covers various positions
        Chem.MolFromSmarts("c1cc(=O)c(=O)c(c1)c2ccccc2"),  # Naphthoquinone-like
        Chem.MolFromSmarts("c1cc2c(c1)c(=O)c3ccccc3c(=O)c2"),  # Anthraquinone-like
        Chem.MolFromSmarts("c1c(=O)c2ccccc2c(=O)c1")  # Covers variations of large conjugations
    ]

    # Check for quinone patterns
    for pattern in quinone_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains fully conjugated cyclic dione structure typical of quinones"
    
    return False, "Does not contain the structural features typical of quinones"