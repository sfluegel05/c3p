"""
Classifies: CHEBI:29256 thiol
"""
"""
Classifies: CHEBI:26806 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol contains at least one -SH group attached to a carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for thiol group (-SH)
    # Matches sulfur atom with exactly one hydrogen (explicit or implicit)
    # attached to any carbon (aliphatic or aromatic)
    thiol_pattern = Chem.MolFromSmarts("[C][SH1]")
    
    # Check for matches
    matches = mol.GetSubstructMatches(thiol_pattern)
    if len(matches) > 0:
        return True, f"Contains {len(matches)} thiol group(s)"
    else:
        return False, "No thiol (-SH) groups detected"