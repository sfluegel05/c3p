"""
Classifies: CHEBI:39362 mononitrophenol
"""
from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is a phenol with a single nitro group attached directly on the same benzene ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mononitrophenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define phenol with a nitro group pattern on the benzene ring
    mononitrophenol_pattern = Chem.MolFromSmarts("c1(ccccc1O)[N+](=O)[O-]")
    
    # Look for a structure that matches exactly one of these mononitrophenol patterns
    matches = mol.GetSubstructMatches(mononitrophenol_pattern)
    if len(matches) == 0:
        return False, "No mononitrophenol structure found"

    # Check that there's exactly one connected phenol and nitro group
    if len(matches) != 1:
        return False, f"Multiple possible mononitrophenol structures found: {len(matches)}"

    return True, "Contains phenol group with a single nitro group on the same benzene ring"