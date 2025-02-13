"""
Classifies: CHEBI:39362 mononitrophenol
"""
from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol must have a benzene ring with at least one hydroxyl group and one nitro group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a mononitrophenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for hydroxyl group anywhere on a benzene ring
    hydroxyl_pattern = Chem.MolFromSmarts("Oc1ccccc1")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No phenol group (hydroxyl on benzene) found"

    # Look for at least one nitro group
    nitro_group_pattern = Chem.MolFromSmarts("[NX3](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_group_pattern)
    if len(nitro_matches) < 1:
        return False, "No nitro group found, need at least 1"

    return True, "Contains at least one hydroxyl group and one nitro group on a benzene ring"