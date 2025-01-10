"""
Classifies: CHEBI:72600 spiroketal
"""
"""
Classifies: spiroketal
A cyclic ketal in which the ketal carbon is the only common atom of two rings.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule contains a spiroketal structure.
    A spiroketal has two rings sharing only one carbon atom (spiro carbon) that is bonded to two oxygens.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains spiroketal, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First pass: look for spiro carbon connected to two oxygens
    # [#6D4] - carbon with 4 bonds
    # (@[#8D2R]) - oxygen in ring with 2 bonds
    # {2} - exactly two occurrences
    spiro_pattern = Chem.MolFromSmarts('[#6D4](@[#8D2R]){2}')
    
    if spiro_pattern is None:
        return False, "Invalid SMARTS pattern"
    
    spiro_matches = mol.GetSubstructMatches(spiro_pattern)
    if not spiro_matches:
        return False, "No carbon atom with two ring oxygens found"

    # For each potential spiro carbon
    for match in spiro_matches:
        spiro_carbon_idx = match[0]
        spiro_carbon = mol.GetAtomWithIdx(spiro_carbon_idx)
        
        # Get ring information
        ring_info = mol.GetRingInfo()
        if not ring_info.NumAtomRings(spiro_carbon_idx):
            continue
            
        # Get all rings containing the spiro carbon
        rings = ring_info.AtomRings()
        rings_with_spiro = []
        
        for ring in rings:
            if spiro_carbon_idx in ring:
                rings_with_spiro.append(set(ring))
        
        # Need at least two rings containing the spiro carbon
        if len(rings_with_spiro) < 2:
            continue
            
        # Check pairs of rings
        for i in range(len(rings_with_spiro)):
            for j in range(i+1, len(rings_with_spiro)):
                ring1 = rings_with_spiro[i]
                ring2 = rings_with_spiro[j]
                
                # The rings should share exactly one atom (the spiro carbon)
                intersection = ring1.intersection(ring2)
                if len(intersection) == 1 and spiro_carbon_idx in intersection:
                    # Verify both rings contain an oxygen
                    ring1_has_o = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring1)
                    ring2_has_o = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring2)
                    
                    if ring1_has_o and ring2_has_o:
                        return True, "Contains spiroketal structure with two rings sharing only one carbon atom bonded to two oxygens"
    
    return False, "No valid spiroketal structure found"