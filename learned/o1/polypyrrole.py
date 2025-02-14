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
    A polypyrrole is a compound composed of two or more pyrrole units.

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

    # Ensure proper aromaticity perception
    Chem.SanitizeMol(mol)

    # Define pyrrole SMARTS pattern
    # Matches a five-membered ring with one nitrogen atom with a hydrogen ([#7&H1])
    # and four carbon atoms, allowing for aromaticity variation and fused rings
    pyrrole_smarts = '[nH]1cccc1'  # Aromatic pyrrole
    pyrrole_pattern = Chem.MolFromSmarts(pyrrole_smarts)

    # Find all matches of the pyrrole substructure
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    pyrrole_count = len(pyrrole_matches)

    # Check for non-aromatic pyrrole units (tautomeric forms)
    if pyrrole_count < 2:
        # Non-aromatic pyrrole pattern
        pyrrole_non_aromatic_smarts = '[#7&H1]-[#6]-[#6]-[#6]-[#6]'
        pyrrole_non_aromatic_pattern = Chem.MolFromSmarts(pyrrole_non_aromatic_smarts)
        pyrrole_non_aromatic_matches = mol.GetSubstructMatches(pyrrole_non_aromatic_pattern)
        pyrrole_count += len(pyrrole_non_aromatic_matches)

    if pyrrole_count >= 2:
        return True, f"Contains {pyrrole_count} pyrrole units"
    else:
        return False, f"Contains only {pyrrole_count} pyrrole unit(s), less than 2 required for polypyrrole"