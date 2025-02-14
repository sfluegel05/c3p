"""
Classifies: CHEBI:38077 polypyrrole
"""
from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole consists of two or more pyrrole units, often manifesting in
    complex ring structures.

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
    
    # Define pyrrole-like patterns for complex systems (fused or extended conjugated systems)
    # More general pyrrole pattern to cover fused systems
    general_pyrrole_pattern = Chem.MolFromSmarts('[nH1]1cccc1 | [n]1ccccc1') 

    # Find matches for pyrrole units
    pyrrole_matches = mol.GetSubstructMatches(general_pyrrole_pattern)
    unique_matches = {tuple(sorted(match)) for match in pyrrole_matches}
    
    # Check if there are two or more distinct pyrrole units
    if len(unique_matches) >= 2:
        return True, f"Contains {len(unique_matches)} pyrrole units"

    return False, "Less than 2 pyrrole units found"