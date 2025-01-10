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

    # 2,5-diketopiperazine pattern:
    # - Six-membered ring
    # - Two nitrogens at positions 1 and 4
    # - Two carbonyls at positions 2 and 5
    dkp_pattern = Chem.MolFromSmarts('[NX3R2]1[CX3R2](=O)[CX4R2][NX3R2][CX3R2](=O)[CX4R2]1')
    
    if dkp_pattern is None:
        return None, "Error in SMARTS pattern"
        
    if not mol.HasSubstructMatch(dkp_pattern):
        return False, "No 2,5-diketopiperazine core structure found"
    
    # Get all matches of the pattern
    matches = mol.GetSubstructMatches(dkp_pattern)
    if len(matches) < 1:
        return False, "No complete 2,5-diketopiperazine ring system found"

    # For each match, verify the structural requirements
    for match in matches:
        # Get the atoms in the match
        ring_atoms = list(match)
        
        # Check that all atoms are in the same ring
        ring_size = 0
        for ring in mol.GetRingInfo().AtomRings():
            if all(idx in ring for idx in ring_atoms):
                ring_size = len(ring)
                break
                
        if ring_size != 6:
            continue  # Not a valid 6-membered ring
            
        # Count the number of carbonyls
        carbonyl_pattern = Chem.MolFromSmarts('[NX3R2][CX3R2](=[OX1])')
        carbonyl_matches = len(mol.GetSubstructMatches(carbonyl_pattern))
        if carbonyl_matches < 2:
            continue  # Not enough carbonyls
            
        # Verify nitrogen atoms are properly placed
        n_atoms = [mol.GetAtomWithIdx(idx) for idx in ring_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7]
        if len(n_atoms) != 2:
            continue  # Wrong number of nitrogens
            
        # All checks passed for this match
        return True, "Contains piperazine-2,5-dione core structure with proper ketone groups"
    
    return False, "No valid 2,5-diketopiperazine structure found"