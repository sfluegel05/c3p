"""
Classifies: CHEBI:38077 polypyrrole
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is a compound composed of two or more pyrrole units.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) indicating if molecule is a polypyrrole and reason
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define multiple pyrrole patterns to catch different cases
    patterns = [
        # Basic pyrrole (both aromatic and non-aromatic)
        '[nH,n]1[c,C][c,C][c,C][c,C]1',
        # Metal-coordinated pyrrole in porphyrins
        '[n]1[c,C][c,C][c,C][c,C]1',
        # Reduced pyrrole (as in bile pigments)
        '[NH,N]1[CH2][CH2][CH2][CH2]1',
        # BODIPY-type structures
        '[n]1[c,C][c,C][c,C][c,C]1[B,Fe,Mg,Zn,Co,Pd]',
        # Fused pyrrole systems
        '[nH,n]1[c,C]2[c,C][c,C][c,C]1[c,C]2',
    ]
    
    total_matches = set()
    for pattern in patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt:
            matches = mol.GetSubstructMatches(patt)
            total_matches.update(matches)
    
    num_pyrroles = len(total_matches)
    
    if num_pyrroles == 0:
        return False, "No pyrrole rings found"
    elif num_pyrroles == 1:
        return False, "Only one pyrrole ring found, need at least two"
    
    # Check for metals that commonly coordinate to pyrroles
    metal_pattern = Chem.MolFromSmarts('[Fe,Mg,Zn,Co,Pd,B]')
    has_metal = mol.HasSubstructMatch(metal_pattern)
    
    # Check for common polypyrrole classes
    porphyrin_core = Chem.MolFromSmarts('[n]1[c,C][c,C][c,C][c,C]1-[c,C]1[c,C][c,C][n]2')
    is_porphyrin_like = mol.HasSubstructMatch(porphyrin_core)
    
    # Construct detailed reason
    reason = f"Contains {num_pyrroles} pyrrole rings"
    if is_porphyrin_like:
        reason += " in a porphyrin-like arrangement"
    if has_metal:
        reason += " with metal coordination"
        
    return True, reason