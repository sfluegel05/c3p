"""
Classifies: CHEBI:27325 xanthophyll
"""
from rdkit import Chem

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    A xanthophyll is characterized by a carotenoid backbone with oxygen functionalities.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved identification of polyene chain for carotenoid backbone
    polyene_pattern = Chem.MolFromSmarts("C=C(-C=C)*C=C")  # More complex and specific polyene pattern
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, "No suitable polyene chain found, not a carotenoid backbone"

    # Look for oxygen atom presence
    has_oxygen = any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())
    if not has_oxygen:
        return False, "No oxygen atoms found, not a xanthophyll"

    # Check for diverse oxygen functionalities relevant to xanthophylls
    oxo_patterns = [
        Chem.MolFromSmarts("[CX3]=[OX1]"),               # Carbonyl
        Chem.MolFromSmarts("[OX2H]"),                    # Hydroxyl
        Chem.MolFromSmarts("[OX2][CX3](=[OX1])[CX4]"),   # Ester-like
        Chem.MolFromSmarts("[OX2][CX3][CX3H1]")          # Epoxide or similar
    ]

    has_oxo_functionality = False
    for pattern in oxo_patterns:
        if mol.HasSubstructMatch(pattern):
            has_oxo_functionality = True
            break

    if not has_oxo_functionality:
        return False, "Lacks key oxygen functionalities like hydroxyl, carbonyl, or other oxygens"

    # All criteria matched, classify as xanthophyll
    return True, "Molecule has a carotenoid backbone and distinctive oxygen functionalities, hence xanthophyll"