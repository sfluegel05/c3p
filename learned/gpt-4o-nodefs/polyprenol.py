"""
Classifies: CHEBI:26199 polyprenol
"""
from rdkit import Chem

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    Polyprenols consist of three or more isoprene units, usually terminating with an alcohol group.

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

    # SMARTS patterns for various forms of isoprene units
    isoprene_patterns = [
        "[C;!R]=[C;!R]-[C;!R]-[C;!R]",  # generic acyclic isoprene-like linking
        "[C;!R]=[C;!R]-[CX4]-[C;!R]",  # with saturated carbon within
        "[CX4]-[C;!R]=[C;!R]-[CX4]"  # fully saturated alternative matching
    ]

    num_isoprene_units = sum(len(mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))) for pattern in isoprene_patterns)
    
    # Polyprenols should have at least 3 isoprene units
    if num_isoprene_units < 3:
        return False, f"Only {num_isoprene_units} isoprene units found, at least 3 required"

    # Check for terminal or near-terminal alcohol groups
    alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H1]")  # alcohol on a saturated carbon
    terminal_alcohol = mol.GetNumAtoms() - 1 in [match[-1] for match in mol.GetSubstructMatches(alcohol_pattern)]

    if not terminal_alcohol:
        return False, "Molecule does not have a terminal alcohol group"

    return True, "Molecule is a polyprenol with three or more isoprene units and a terminal alcohol group"