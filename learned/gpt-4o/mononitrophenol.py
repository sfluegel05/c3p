"""
Classifies: CHEBI:39362 mononitrophenol
"""
from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is a phenol (benzene with an -OH group) carrying a single nitro substituent attached to the aromatic ring.
    
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
    
    # Look for phenol group (benzene ring with -OH)
    phenol_pattern = Chem.MolFromSmarts("c1ccccc1O")
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    if len(phenol_matches) == 0:
        return False, "No phenol group found"
    
    # Look for nitro groups attached to aromatic carbon
    nitro_pattern = Chem.MolFromSmarts("c[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    
    # Check nitro group count
    if len(nitro_matches) != 1:
        return False, f"Expected 1 nitro group attached, found {len(nitro_matches)}"
    
    # Ensure nitro is attached to the phenol's aromatic ring
    # Convert matches to sets of atom indices for easier comparison
    phenol_indices = set(sum(phenol_matches, ()))
    nitro_indices = set(nitro_matches[0])
    
    # Check intersection: there must be an overlap in indices
    if len(phenol_indices.intersection(nitro_indices)) == 0:
        return False, "Nitro group is not attached to the phenol ring"

    return True, "Contains phenol group with exactly one nitro group attached"