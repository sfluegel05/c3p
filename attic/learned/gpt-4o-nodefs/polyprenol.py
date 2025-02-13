"""
Classifies: CHEBI:26199 polyprenol
"""
from rdkit import Chem

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    Polyprenols are oligomers with three or more isoprene units, often ending in
    an alcohol group.

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

    # SMARTS pattern for generic isoprene unit in various configurations: CH=C-C=C or C=C-C=C
    # This accommodates variations in isoprene subunit connectivity seen in polyprenols
    isoprene_smarts = "[C;R0]=[C;R0]-[C;R0]-[C;R0]"

    isoprene_pattern = Chem.MolFromSmarts(isoprene_smarts)

    # Find substructure matches
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)

    # Polyprenols should have at least 3 such isoprene units
    if len(isoprene_matches) < 3:
        return False, f"Only {len(isoprene_matches)} isoprene units found, at least 3 required"

    # Check for alcohol group (OH) at the terminal
    # This pattern looks for terminal or near-terminal OH bound to a saturated carbon (which is often Csp3 in alcohols)
    alcohol_pattern = Chem.MolFromSmarts("[OX2H1]")  # Hydroxyl group
    
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "Molecule does not have a terminal alcohol group"

    return True, "Molecule is a polyprenol with three or more isoprene units and a terminal alcohol group"