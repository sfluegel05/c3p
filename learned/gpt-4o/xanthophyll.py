"""
Classifies: CHEBI:27325 xanthophyll
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Identify if there is a polyene chain (long conjugated chain)
    polyene_pattern = Chem.MolFromSmarts("C=C[C]=C") # Slightly broader pattern
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, "No polyene chain found, not a carotenoid backbone"

    # Look for oxygenated functionalities
    has_oxygen = any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())
    if not has_oxygen:
        return False, "No oxygen atoms found, not a xanthophyll"

    # Check for specific oxygen functionalities
    oxo_patterns = [
        Chem.MolFromSmarts("[CX3]=[OX1]"),  # Carbonyl group
        Chem.MolFromSmarts("[OX2H]"),       # Hydroxyl group
        Chem.MolFromSmarts("[OX2]([CX3]=[OX1])"),  # Epoxide groups
        Chem.MolFromSmarts("[O][CX3]([O])=O")      # Additional functional forms
    ]

    has_oxo_functionality = any(mol.HasSubstructMatch(pattern) for pattern in oxo_patterns)
    if not has_oxo_functionality:
        return False, "No oxygen functionalities like hydroxyl, carbonyl or epoxide groups found"

    return True, "Molecule has a carotenoid backbone and oxygen functionalities, classifying it as xanthophyll"