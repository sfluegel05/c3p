"""
Classifies: CHEBI:38077 polypyrrole
"""
from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    Polypyrroles generally consist of pyrrole rings, but are ambiguously defined here.

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
    
    # Look for pyrrole-like substructure
    # Pyrrole ring is C1=CNC=C1, using SMARTS pattern for pyrrole
    pyrrole_pattern = Chem.MolFromSmarts('n1cccc1')
    if mol.HasSubstructMatch(pyrrole_pattern):
        # Since specific definition is not given, assume any molecule containing pyrrole ring is polypyrrole
        return True, "Contains pyrrole ring"

    return False, "No pyrrole structure found or definition for polypyrrole is 'None'"