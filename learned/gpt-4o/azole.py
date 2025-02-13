"""
Classifies: CHEBI:68452 azole
"""
from rdkit import Chem

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole based on its SMILES string.
    An azole is defined as a monocyclic or fused five-membered ring (heteroarene) containing nitrogen and possible other heteroatoms such as N, S, or O.

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

    # Extended SMARTS patterns to account for different azole forms
    azole_patterns = [
        Chem.MolFromSmarts("n1cccn1"),  # Imidazole
        Chem.MolFromSmarts("n1c[nH]cc1"),  # Pyrazole
        Chem.MolFromSmarts("n1ncc[nH]1"),  # 1,2,3-Triazole
        Chem.MolFromSmarts("n1nncn1"),  # Tetrazole
        Chem.MolFromSmarts("c1ncn[nH]1"),  # 1,2,4-Triazole
        Chem.MolFromSmarts("n1cn[o,s]c1"),  # Oxazole/Thiazole
        Chem.MolFromSmarts("n1csc[nH]1"),  # Thiazole variant
        Chem.MolFromSmarts("n1onc2ccccc2c1"), # Furan-like plus azole in benzo-fused systems
        Chem.MolFromSmarts("c1nccn1"),  # Pyrrole-like azole
    ]

    # Check if any azole pattern is present in the molecule
    for pattern in azole_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains an azole ring based on expanded set of patterns"

    return False, "No azole ring found"