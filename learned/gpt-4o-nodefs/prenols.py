"""
Classifies: CHEBI:26244 prenols
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    Prenols generally have a terminal alcohol and consist of multiple connected isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a terminal alcohol group '-OH' at the end of the carbon chain
    # Ensure it does not match anywhere else arbitrarily
    terminal_oh_pattern = Chem.MolFromSmarts("[CX4,!R]-[OX2H]")
    terminal_oh_matches = mol.GetSubstructMatches(terminal_oh_pattern)
    if len(terminal_oh_matches) == 0:
        return False, "Missing terminal alcohol group"

    # Look for sequences of isoprene units: characterized by multiple C(=C-C-C=C) patterns
    # Use recursive SMARTS pattern to match sequences
    isoprene_pattern = Chem.MolFromSmarts("(C(=C)C-C(=C)C)")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 1:
        return False, "No significant isoprene unit sequences detected"
    
    return True, "Contains terminal alcohol group and isoprene units, typical of prenols"