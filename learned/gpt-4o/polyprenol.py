"""
Classifies: CHEBI:26199 polyprenol
"""
from rdkit import Chem

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol has more than one isoprene unit and ends with an alcohol group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for isoprene unit: C=C-C-C
    isoprene_pattern = Chem.MolFromSmarts("C(=C)C-C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    
    # Check if there is more than one isoprene unit
    if len(isoprene_matches) < 2:
        return False, f"Found {len(isoprene_matches)} isoprene units, need more than one"

    # Pattern for terminal alcohol group
    alcohol_pattern = Chem.MolFromSmarts("[OX2H1]")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No terminal alcohol group found"

    return True, "Contains more than one isoprene unit and a terminal alcohol group"