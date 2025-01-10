"""
Classifies: CHEBI:38077 polypyrrole
"""
"""
Classifies: polypyrrole
"""
from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is defined as a compound composed of two or more pyrrole units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypyrrole, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a broader pyrrole pattern to include variations
    pyrrole_patterns = [
        Chem.MolFromSmarts("c1[nH]ccc1"),  # Common pyrrole
        Chem.MolFromSmarts("n1cccc1"),    # Pyrrole with different bond saturation
        Chem.MolFromSmarts("[#7]1[#6]=[#6][#6]=[#6]1")  # Capturing possible tautomers/form variants
    ]
    
    # Collect all matches from all patterns
    pyrrole_matches = set()
    for pattern in pyrrole_patterns:
        pyrrole_matches.update(mol.GetSubstructMatches(pattern))
    
    num_pyrrole_units = len(pyrrole_matches)
    
    if num_pyrrole_units < 2:
        return False, f"Contains only {num_pyrrole_units} pyrrole units, need 2 or more for polypyrrole"

    return True, f"Contains {num_pyrrole_units} pyrrole units, qualifies as polypyrrole"