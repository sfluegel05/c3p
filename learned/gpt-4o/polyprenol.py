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

    # Define isoprene unit pattern: CH2=C-CH=CH-CH2
    isoprene_pattern = Chem.MolFromSmarts("C(=C)C=C-C")
    
    # Find isoprene units in the molecule
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 2:
        return False, f"Found {len(isoprene_matches)} isoprene units, need more than one"

    # Ensure terminal alcohol group - simplified pattern
    terminal_alcohol_pattern = Chem.MolFromSmarts("[CX4][OH]")
    alcohol_matches = mol.GetSubstructMatches(terminal_alcohol_pattern)
    if not any(mol.GetAtomWithIdx(match[1]).GetSymbol() == 'O' for match in alcohol_matches):
        return False, "No terminal alcohol group found"
    
    return True, "Contains more than one isoprene unit and a terminal alcohol group"