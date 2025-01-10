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

    # More flexible 2,5-diketopiperazine pattern:
    # - Six-membered ring with any bond types
    # - Two nitrogens
    # - Two carbonyls
    # Note: removed specific hybridization requirements and bond restrictions
    dkp_pattern = Chem.MolFromSmarts('[N]1[C](=O)[C][N][C](=O)[C]1')
    
    if dkp_pattern is None:
        return None, "Error in SMARTS pattern"
        
    # Alternative pattern for cases with double bonds in ring
    dkp_pattern2 = Chem.MolFromSmarts('[N]1[C](=O)[C]=[N][C](=O)[C]1')
    
    # Pattern for bridged/fused systems
    dkp_pattern3 = Chem.MolFromSmarts('[N]1[C](=O)[C]2[N][C](=O)[C]1[*]2')
    
    if not (mol.HasSubstructMatch(dkp_pattern) or 
            mol.HasSubstructMatch(dkp_pattern2) or 
            mol.HasSubstructMatch(dkp_pattern3)):
        return False, "No 2,5-diketopiperazine core structure found"
    
    # Verify the presence of two carbonyls adjacent to nitrogens
    carbonyl_pattern = Chem.MolFromSmarts('[N][C](=O)')
    carbonyl_matches = len(mol.GetSubstructMatches(carbonyl_pattern))
    
    if carbonyl_matches < 2:
        return False, "Missing required carbonyl groups"
        
    # Additional check for ring system
    ring_found = False
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) >= 6:  # Allow for fused systems
            # Check if ring contains required elements
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            n_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 7)
            if n_count >= 2:
                ring_found = True
                break
    
    if not ring_found:
        return False, "No suitable ring system found"

    return True, "Contains piperazine-2,5-dione core structure"