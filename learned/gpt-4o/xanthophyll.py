"""
Classifies: CHEBI:27325 xanthophyll
"""
from rdkit import Chem

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    A xanthophyll is characterized by a conjugated polyene carbons backbone (typical of carotenoids)
    and is oxygenated.

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

    # Identify conjugated polyene chain pattern
    polyene_chain_pattern = Chem.MolFromSmarts("C=CC=CC=C")  # Basic pattern for conjugated dienes
    if not mol.HasSubstructMatch(polyene_chain_pattern):
        return False, "No suitable conjugated polyene chain found, not a carotenoid backbone"
    
    # Xanthophyll specific oxygen-containing groups presence check
    oxygen_pattern = Chem.MolFromSmarts("[OX2H1]")  # Presence of hydroxyl group
    if not mol.HasSubstructMatch(oxygen_pattern):
        return False, "No hydroxyl functionalities found, possibly missing oxygen functionalities"

    # Look for oxygen atoms indicating potential other functionalities
    has_oxygen = any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())
    if not has_oxygen:
        return False, "No oxygen atoms found, not a xanthophyll"

    # Additional check: specific xanthophyll oxygen functionalities (e.g., epoxide)
    oxo_patterns = [
        Chem.MolFromSmarts("[CX3]=[OX1]"),               # Carbonyl indicates possible ketone
        Chem.MolFromSmarts("[OX2][CX3](=[OX1])"),        # Ester-like (common in modifications)
        Chem.MolFromSmarts("[CX4][OX2][CX4]")            # Epoxide group presence check
    ]
    
    has_oxo_functionality = any(mol.HasSubstructMatch(pattern) for pattern in oxo_patterns)
    if not has_oxo_functionality:
        return False, "Lacks key oxygen functionalities like hydroxyl or carbonyl"

    # If all criteria matched, classify as xanthophyll
    return True, "Molecule has a carotenoid backbone with conjugated polyene and oxygen functionalities, hence xanthophyll"