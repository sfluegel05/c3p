"""
Classifies: CHEBI:26199 polyprenol
"""
from rdkit import Chem

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol consists of more than one isoprene unit and ends with an alcohol group.

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

    # Define a more flexible isoprene unit pattern
    # Covers variations including potential branching and steric configurations
    isoprene_pattern = Chem.MolFromSmarts("C(=C)CC=C[CH2]") # General pattern allowing variations in branching

    # Verify isoprene unit matches
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 2:
        return False, f"Found {len(isoprene_matches)} isoprene units, need more than one"

    # Define a generalized pattern for terminal alcohol group
    # Detect primary alcohols at the end of the carbon chain
    terminal_alcohol = Chem.MolFromSmarts("[CX4;R0][OH]")

    # Verify terminal alcohol matches
    alcohol_matches = mol.GetSubstructMatches(terminal_alcohol)
    if not any(mol.GetAtomWithIdx(match[1]).GetSymbol() == 'O' for match in alcohol_matches):
        return False, "No terminal alcohol group found"
    
    return True, "Contains more than one isoprene unit and a terminal alcohol group"