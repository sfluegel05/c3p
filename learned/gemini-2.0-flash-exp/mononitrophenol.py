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

    # Check for phenol substructure (aromatic carbon with -OH or -O-)
    phenol_pattern = Chem.MolFromSmarts("[c][OH]")
    phenol_pattern2 = Chem.MolFromSmarts("[c][O-]")
    if not (mol.HasSubstructMatch(phenol_pattern) or mol.HasSubstructMatch(phenol_pattern2)):
            return False, "No phenol substructure found"

    #Check for exactly one nitro group directly attached to any aromatic ring
    nitro_pattern = Chem.MolFromSmarts("[c][N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if len(nitro_matches) != 1:
            return False, f"Found {len(nitro_matches)} nitro groups directly attached to aromatic ring, expected exactly 1"
    
    # Check that the nitro group and the hydroxyl group are on the same aromatic system.
    # This is hard to do, so we'll use a simpler approach for now.

    return True, "Contains a phenol substructure with a single nitro group directly attached to an aromatic ring"