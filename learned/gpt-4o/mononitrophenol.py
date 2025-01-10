"""
Classifies: CHEBI:39362 mononitrophenol
"""
from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is a phenol carrying a single nitro substituent.

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

    # Look for phenol group pattern
    phenol_pattern = Chem.MolFromSmarts("c1ccc(cc1)O")
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    if len(phenol_matches) != 1:
        return False, "Phenol group not found or multiple phenol groups present"
        
    # Look for nitro group pattern
    nitro_pattern = Chem.MolFromSmarts("[N+]([O-])=O")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if len(nitro_matches) != 1:
        return False, f"Found {len(nitro_matches)} nitro groups, need exactly 1"

    return True, "Contains single phenol group with a single nitro group attached"