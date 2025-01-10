"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
from rdkit import Chem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    A hydroxynaphthoquinone is a naphthoquinone moiety with at least one hydroxy group substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxynaphthoquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for naphthoquinone: aromatic naphthalene with two ketones
    naphthoquinone_pattern = Chem.MolFromSmarts('O=C1C=Cc2ccccc2C1=O')
    if naphthoquinone_pattern is None:
        return False, "Invalid naphthoquinone SMARTS pattern"

    # SMARTS pattern for a hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts('[OH]')
    if hydroxy_pattern is None:
        return False, "Invalid hydroxy SMARTS pattern"

    # Check for naphthoquinone core structure
    if not mol.HasSubstructMatch(naphthoquinone_pattern):
        return False, "No naphthoquinone core found"

    # Check for hydroxy group presence
    if mol.HasSubstructMatch(hydroxy_pattern):
        return True, "Contains naphthoquinone core with hydroxy group substitution"

    return False, "Naphthoquinone core present but no hydroxy substitution found"