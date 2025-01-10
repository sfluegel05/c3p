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
    
    # Define phenol pattern (benzene ring with -OH group)
    phenol_pattern = Chem.MolFromSmarts("c1ccccc1O")
    # Define nitro group pattern
    nitro_pattern = Chem.MolFromSmarts("[$(a-[N+](=O)[O-])]")
    
    # Check if phenol structure is present
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    if len(phenol_matches) == 0:
        return False, "No phenol group found"
    
    # Check for nitro group
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if len(nitro_matches) != 1:
        return False, f"Found {len(nitro_matches)} nitro groups, need exactly 1"
    
    # Verify nitro is directly bonded to the aromatic carbon of phenol
    for match in phenol_matches:
        benzene_atoms = set(match[:-1])  # exclude the OH
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if ({begin_idx, end_idx} & benzene_atoms) and (begin_idx in nitro_matches[0] or end_idx in nitro_matches[0]):
                return True, "Molecule is a mononitrophenol with a nitro group attached to the phenol ring"
    
    return False, "Nitro group is not attached to the aromatic ring of phenol"