"""
Classifies: CHEBI:39362 mononitrophenol
"""
import itertools
from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol must have at least one hydroxyl group and one nitro group attached to the same aromatic system.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a mononitrophenol, False otherwise.
        str: Reason for classification.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Check for phenolic hydroxyl in an aromatic ring
    hydroxyl_pattern = Chem.MolFromSmarts("c1cc(ccc1)O")  # This ensures hydroxyl is part of an aromatic system
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No phenolic hydroxyl group found on aromatic system"
    
    # Find nitro group pattern
    nitro_group_pattern = Chem.MolFromSmarts("[NX3](=O)[O-]")  # Nitro group pattern
    nitro_matches = mol.GetSubstructMatches(nitro_group_pattern)
    if not nitro_matches:
        return False, "No nitro group found, need at least one nitro group"
    
    # Analyze if nitro and hydroxyl groups are on the same aromatic ring
    hydroxyl_atoms = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_atoms:
        return False, "No phenolic hydroxyl group found on aromatic system"
    
    # Get ring information
    ri = mol.GetRingInfo()
    for match in itertools.product(hydroxyl_atoms, nitro_matches):
        hydroxyl_idx = match[0][0]
        nitro_idx = match[1][0]
        
        # Check if both are in the same ring
        for ring in ri.AtomRings():
            if hydroxyl_idx in ring and nitro_idx in ring:
                return True, "Phenolic hydroxyl and nitro group on the same aromatic system found"
    
    return False, "Phenolic hydroxyl and nitro group not on the same benzene ring"