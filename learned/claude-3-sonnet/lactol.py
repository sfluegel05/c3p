"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies: CHEBI:35879 lactol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    Lactols are cyclic hemiacetals formed by intramolecular addition of a hydroxy 
    group to an aldehydic or ketonic carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cyclic hemiacetal pattern: O-C-OH where O is in ring
    # [O;R] means oxygen in ring
    # [C;R] means carbon in ring
    # [OH1] means hydroxyl group
    lactol_pattern = Chem.MolFromSmarts("[O;R][C;R][OH1]")
    
    # Alternative pattern for lactol where the OH is in ring
    lactol_pattern2 = Chem.MolFromSmarts("[O;R][C;R][O;R]")
    
    matches = mol.GetSubstructMatches(lactol_pattern)
    matches2 = mol.GetSubstructMatches(lactol_pattern2)
    
    if not (matches or matches2):
        return False, "No cyclic hemiacetal (O-C-OH) structure found"

    # Additional check to filter out false positives
    # Count number of unique rings that contain the lactol pattern
    ring_info = mol.GetRingInfo()
    
    valid_lactol = False
    for match in matches + matches2:
        # Check if oxygen and carbon are in same ring
        o_atom, c_atom = match[0], match[1]
        rings_with_o = set(ring_info.AtomRings()[i] 
                          for i, ring in enumerate(ring_info.AtomRings()) 
                          if o_atom in ring)
        rings_with_c = set(ring_info.AtomRings()[i] 
                          for i, ring in enumerate(ring_info.AtomRings()) 
                          if c_atom in ring)
        
        # If O and C share at least one ring and pattern matches,
        # we have found a valid lactol
        if rings_with_o.intersection(rings_with_c):
            valid_lactol = True
            break

    if not valid_lactol:
        return False, "Hemiacetal oxygen and carbon must be in same ring"

    return True, "Contains cyclic hemiacetal structure (lactol)"