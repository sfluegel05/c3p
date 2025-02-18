"""
Classifies: CHEBI:83970 cardiac glycoside
"""
"""
Classifies: CHEBI:37701 cardiac glycoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    Cardiac glycosides contain a steroid nucleus with a lactone ring and sugar residues.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a cardiac glycoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # 1. Check for butenolide lactone ring (O=C1OC=CC1)
    lactone_pattern = Chem.MolFromSmarts('O=C1OC=CC1')
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    if not lactone_matches:
        return False, "No lactone ring detected"
    
    # 2. Verify lactone is attached to fused ring system (steroid nucleus)
    fused_ring_found = False
    for match in lactone_matches:
        # Get carbon attached to lactone oxygen (position 1 in SMARTS)
        c_idx = match[1]
        c_atom = mol.GetAtomWithIdx(c_idx)
        # Check if this carbon is in at least two fused rings
        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()
        ring_count = sum(1 for ring in rings if c_idx in ring)
        if ring_count >= 2:
            fused_ring_found = True
            break
    if not fused_ring_found:
        return False, "Lactone not attached to fused ring system (steroid)"
    
    # 3. Check for sugar residues via glycosidic bonds
    # Find oxygen atoms in rings (potential sugar moieties)
    sugar_o_pattern = Chem.MolFromSmarts('[O;R]')
    sugar_matches = mol.GetSubstructMatches(sugar_o_pattern)
    if not sugar_matches:
        return False, "No sugar residues found"
    
    # Verify glycosidic bond connects sugar to steroid/lactone core
    glycosidic_bond = False
    for o_match in sugar_matches:
        o_atom = mol.GetAtomWithIdx(o_match[0])
        # Check neighbors for connection to non-sugar part
        for neighbor in o_atom.GetNeighbors():
            # Check if neighbor atom is not in a sugar ring
            in_sugar = False
            for ring in rings:
                if neighbor.GetIdx() in ring:
                    # Count oxygens in ring to determine sugar likelihood
                    o_count = sum(1 for a in ring if mol.GetAtomWithIdx(a).GetAtomicNum() == 8)
                    if o_count >= 3:  # Heuristic for sugar rings
                        in_sugar = True
                        break
            if not in_sugar:
                glycosidic_bond = True
                break
        if glycosidic_bond:
            break
    if not glycosidic_bond:
        return False, "No glycosidic bond linking sugar to core"
    
    return True, "Steroid lactone with glycosidic sugar residues"