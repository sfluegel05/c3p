"""
Classifies: CHEBI:39362 mononitrophenol
"""
from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is defined as a nitrophenol where the phenol carries a
    single nitro substituent at an unspecified position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mononitrophenol, False otherwise
        tuple: Reason for classification or potential compound details
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define phenol pattern (aromatic ring with -OH group)
    phenol_pattern = Chem.MolFromSmarts("c1ccc(cc1)O")
    # Define nitro group pattern
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    
    # Check if phenol structure is present
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    if len(phenol_matches) == 0:
        return False, "No phenol group found"
    
    # Check for nitro group and its attachment to the aromatic ring
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if len(nitro_matches) != 1:
        return False, f"Found {len(nitro_matches)} nitro groups, need exactly 1"
    
    # Ensure at least one match within aromatic system
    phenol_aromatic_atoms = {atom_idx for match in phenol_matches for atom_idx in match}
    if not any(atom_idx in phenol_aromatic_atoms for atom_idx, _, _ in nitro_matches):
        return False, "Nitro group is not attached to the aromatic ring of phenol"
    
    return True, "Molecule is a mononitrophenol with a single nitro group attached to the phenol ring"