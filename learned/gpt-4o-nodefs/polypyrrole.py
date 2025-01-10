"""
Classifies: CHEBI:38077 polypyrrole
"""
from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is a conductive polymer composed of repeating pyrrole units linked via carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as polypyrrole, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define an improved pattern for polymeric pyrrole chain (linking via carbon atoms)
    # The pattern looks for two pyrrole rings (n1cccc1) linked through any carbon atom.
    pyrrole_unit = 'n1cccc1'
    polypyrrole_pattern = Chem.MolFromSmarts(f'{pyrrole_unit}[$([#6]!@[#6]-{pyrrole_unit})]')
    
    # Exclude macrocycles and porphyrin-like structures in a more generic way
    cyclic_pyrrole_exclusion = Chem.MolFromSmarts('c1cc[nH]c([cR3])[nH]c1')  # Extended to catch common cyclic patterns

    # Check for the presence of polypyrrole pattern and absence of large cyclic structures
    if mol.HasSubstructMatch(polypyrrole_pattern) and not mol.HasSubstructMatch(cyclic_pyrrole_exclusion):
        return True, "Contains a polymeric chain of pyrrole units connected via carbon atoms"

    # If pyrrole is present but not forming a polypyrrole
    if mol.HasSubstructMatch(Chem.MolFromSmarts(pyrrole_unit)):
        return False, "Contains pyrrole ring, but not in the form of a polypyrrole chain"

    return False, "No polypyrrole-defining structure found: lacks polymeric pyrrole chain"