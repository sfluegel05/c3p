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

    # Look for basic cyclic hemiacetal pattern in 5 or 6-membered rings
    # [O;R{5-6}] means oxygen in 5 or 6-membered ring
    # [C;R{5-6}] means carbon in 5 or 6-membered ring
    # [OH1] means hydroxyl group
    lactol_pattern = Chem.MolFromSmarts("[O;R{5-6}][C;R{5-6}][OH1]")
    
    # Alternative pattern for lactol where both oxygens are in ring
    lactol_pattern2 = Chem.MolFromSmarts("[O;R{5-6}][C;R{5-6}][O;R{5-6}]")
    
    # Pattern for lactol in equilibrium with open chain form
    # (captures cases where the carbonyl is more prominent)
    lactol_pattern3 = Chem.MolFromSmarts("[O;R{5-6}][C;R{5-6}](=O)")
    
    # Pattern to exclude glycosides (sugar-like structures with OR at anomeric carbon)
    glycoside_pattern = Chem.MolFromSmarts("[OR0][C;R]1[C;R][C;R][C;R][C;R][C;R]1")
    
    matches = mol.GetSubstructMatches(lactol_pattern)
    matches2 = mol.GetSubstructMatches(lactol_pattern2)
    matches3 = mol.GetSubstructMatches(lactol_pattern3)
    glycoside_matches = mol.GetSubstructMatches(glycoside_pattern)
    
    if not (matches or matches2 or matches3):
        return False, "No lactol structure found"

    # Get ring information
    ring_info = mol.GetRingInfo()
    
    # Check for valid lactol structure
    valid_lactol = False
    for match in matches + matches2 + matches3:
        # Check if atoms are in same ring
        o_atom, c_atom = match[0], match[1]
        
        # Get rings containing these atoms
        rings_with_o = set(ring_info.AtomRings()[i] 
                          for i, ring in enumerate(ring_info.AtomRings()) 
                          if o_atom in ring)
        rings_with_c = set(ring_info.AtomRings()[i] 
                          for i, ring in enumerate(ring_info.AtomRings()) 
                          if c_atom in ring)
        
        # Find common rings
        common_rings = rings_with_o.intersection(rings_with_c)
        
        # Check ring size (should be 5 or 6)
        for ring in common_rings:
            if len(ring) in [5, 6]:
                valid_lactol = True
                break
                
        if valid_lactol:
            break

    if not valid_lactol:
        return False, "No valid lactol ring structure found"

    # Exclude molecules that appear to be primarily glycosides
    if len(glycoside_matches) > 0 and len(matches + matches2 + matches3) == 1:
        return False, "Structure appears to be a glycoside rather than a lactol"

    return True, "Contains cyclic hemiacetal structure (lactol)"