"""
Classifies: CHEBI:39362 mononitrophenol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is a phenol with a single nitro substituent at any position.

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

    # Check for phenol substructure
    phenol_pattern = Chem.MolFromSmarts("c1ccccc1O")
    if not mol.HasSubstructMatch(phenol_pattern):
        return False, "No phenol substructure found"
    
    # Check for exactly one nitro group
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if len(nitro_matches) != 1:
         return False, f"Found {len(nitro_matches)} nitro groups, expected exactly 1"

    return True, "Contains a phenol substructure with a single nitro group"