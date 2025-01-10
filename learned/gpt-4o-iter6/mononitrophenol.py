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
    phenol_pattern = Chem.MolFromSmarts("aO")
    # Define nitro group pattern
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    
    # Check if phenol structure is present
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    if len(phenol_matches) == 0:
        return False, "No phenol group found"
    
    # Check for nitro group
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if len(nitro_matches) != 1:
        return False, f"Found {len(nitro_matches)} nitro groups, need exactly 1"
    
    # Check if the nitro group is directly bonded to the aromatic ring
    phenol_aromatic_atoms = {atom_idx for match in phenol_matches for atom_idx in match}
    nitro_atom_idx, _, _ = nitro_matches[0]
    
    # Ensure nitro is bonded to any aromatic carbon (part of phenol)
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() == nitro_atom_idx or bond.GetEndAtomIdx() == nitro_atom_idx:
            if {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()} & phenol_aromatic_atoms:
                return True, "Molecule is a mononitrophenol with a single nitro group attached to the phenol ring"
    
    return False, "Nitro group is not attached to the aromatic ring of phenol"