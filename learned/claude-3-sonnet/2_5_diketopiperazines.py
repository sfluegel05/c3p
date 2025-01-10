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

    # Basic 2,5-diketopiperazine pattern with explicit hydrogens removed
    # Matches the core six-membered ring with two nitrogens and two ketone groups
    dkp_pattern = Chem.MolFromSmarts("[NR2]-[CR2](=O)-[CR2]-[NR2]-[CR2](=O)-[CR2]")
    
    if dkp_pattern is None:
        return None, "Error in SMARTS pattern"
        
    if not mol.HasSubstructMatch(dkp_pattern):
        return False, "No 2,5-diketopiperazine core structure found"
    
    # Count the number of matches to ensure we have the right pattern
    matches = mol.GetSubstructMatches(dkp_pattern)
    if len(matches) < 1:
        return False, "No complete 2,5-diketopiperazine ring system found"
        
    # Verify the presence of two ketone groups connected to nitrogens
    ketone_pattern = Chem.MolFromSmarts("[NR2]-[CR2](=[OX1])")
    if ketone_pattern is None:
        return None, "Error in ketone SMARTS pattern"
        
    ketone_matches = len(mol.GetSubstructMatches(ketone_pattern))
    if ketone_matches < 2:
        return False, f"Found only {ketone_matches} ketone groups, need exactly 2"

    # Additional check for the ring system
    ring_info = mol.GetRingInfo()
    if not any(size == 6 for size in ring_info.RingSizes()):
        return False, "No six-membered ring found"

    # Verify that the nitrogens are part of the same ring
    ring_atoms = mol.GetSubstructMatch(dkp_pattern)
    if not ring_atoms:
        return False, "Ring atoms not properly connected"
        
    # Get the ring bond count for each nitrogen to ensure they're in the ring
    for atom_idx in ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 7:  # Nitrogen
            if not atom.IsInRing():
                return False, "Nitrogen atoms must be part of the ring"

    return True, "Contains piperazine-2,5-dione core structure with proper ketone groups"