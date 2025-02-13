"""
Classifies: CHEBI:39362 mononitrophenol
"""
from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is a phenol with a single nitro group attached at unspecified position on the benzene ring.

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

    # Look for phenol pattern (benzene ring with hydroxyl group)
    phenol_pattern = Chem.MolFromSmarts("c1ccccc1O")
    if not mol.HasSubstructMatch(phenol_pattern):
        return False, "No phenol group found"

    # Look for single nitro group (NO2)
    nitro_pattern = Chem.MolFromSmarts("[N+]([O-])=O")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if len(nitro_matches) != 1:
        return False, f"Found {len(nitro_matches)} nitro groups, need exactly 1"

    return True, "Contains phenol group with a single nitro group at unspecified position"