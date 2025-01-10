"""
Classifies: CHEBI:39362 mononitrophenol
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is a phenol (any benzene ring with an -OH group) carrying a single
    nitro substituent attached to any position on the aromatic ring.

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
    phenol_pattern = Chem.MolFromSmarts("c1c(O)cccc1")  # More flexible aromatic carbon with hydroxyl
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    if len(phenol_matches) == 0:
        return False, "No phenol group found"
    
    # Look for nitro groups
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")  # General search for nitro group
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    
    # Check nitro group count
    if len(nitro_matches) != 1:
        return False, f"Expected 1 nitro group attached, found {len(nitro_matches)}"
    
    # Ensure nitro is attached to the aromatic ring
    for match in phenol_matches:
        aromatic_ring = set(match)
        for nitro in nitro_matches:
            nitro_attachment = set(nitro)
            if any(mol.GetBondBetweenAtoms(i, j) is not None \
                   and mol.GetBondBetweenAtoms(i, j).IsAromatic() for i in aromatic_ring for j in nitro_attachment):
                return True, "Contains phenol group with exactly one nitro group attached"

    return False, "Nitro group is not attached to the phenol ring"