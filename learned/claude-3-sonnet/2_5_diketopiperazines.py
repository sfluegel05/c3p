"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: 2,5-diketopiperazines
Definition: Any piperazinone that has a piperazine-2,5-dione skeleton
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (True if molecule is a 2,5-diketopiperazine, reason for classification)
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core 2,5-diketopiperazine pattern:
    # Two nitrogens connected by carbons, with ketone groups at positions 2 and 5
    dkp_pattern = Chem.MolFromSmarts("[NX3R2]-[CX3R2](=[OX1])-[CX4R2]-[NX3R2]-[CX3R2](=[OX1])-[CX4R2]-1")
    
    # Alternative pattern that can catch some variations
    dkp_pattern2 = Chem.MolFromSmarts("[NR2]-[CR2](=O)-[CR2]-[NR2]-[CR2](=O)-[CR2]")
    
    if not (mol.HasSubstructMatch(dkp_pattern) or mol.HasSubstructMatch(dkp_pattern2)):
        return False, "No 2,5-diketopiperazine core structure found"

    # Check for the presence of exactly two ketone groups in the ring
    ketone_pattern = Chem.MolFromSmarts("[CX3R2](=[OX1])-[NX3R2]")
    ketone_matches = len(mol.GetSubstructMatches(ketone_pattern))
    
    if ketone_matches < 2:
        return False, f"Found only {ketone_matches} ketone groups, need at least 2"

    # Check for six-membered ring with alternating carbons and nitrogens
    ring_pattern = Chem.MolFromSmarts("[NR2]-[CR2]-[CR2]-[NR2]-[CR2]-[CR2]")
    if not mol.HasSubstructMatch(ring_pattern):
        return False, "No proper six-membered ring structure found"

    # Additional check for ring size
    ring_info = mol.GetRingInfo()
    has_six_membered_ring = False
    for ring_size in ring_info.RingSizes():
        if ring_size == 6:
            has_six_membered_ring = True
            break
            
    if not has_six_membered_ring:
        return False, "No six-membered ring found"

    return True, "Contains piperazine-2,5-dione core structure with proper ketone groups"