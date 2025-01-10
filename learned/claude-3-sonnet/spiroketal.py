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

    # Look for potential spiro carbons
    # [#6] is carbon
    # [#8] is oxygen
    # Pattern looks for carbon with two ring oxygens attached
    spiro_pattern = Chem.MolFromSmarts('[#6]([#8;R])[#8;R]')
    
    if not mol.HasSubstructMatch(spiro_pattern):
        return False, "No carbon atom with two ring oxygens found"
    
    matches = mol.GetSubstructMatches(spiro_pattern)
    
    # For each potential spiro carbon
    for match in matches:
        spiro_carbon_idx = match[0]
        o1_idx = match[1]
        o2_idx = match[2]
        
        # Get ring information
        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()
        
        # Find rings containing each oxygen
        rings_o1 = []
        rings_o2 = []
        
        for ring in rings:
            if o1_idx in ring:
                rings_o1.append(set(ring))
            if o2_idx in ring:
                rings_o2.append(set(ring))
        
        # Check each combination of rings
        for ring1 in rings_o1:
            for ring2 in rings_o2:
                # The rings should share exactly one atom (the spiro carbon)
                intersection = ring1.intersection(ring2)
                if len(intersection) == 1 and spiro_carbon_idx in intersection:
                    # Verify the carbon is sp3 (4 bonds)
                    spiro_carbon = mol.GetAtomWithIdx(spiro_carbon_idx)
                    if spiro_carbon.GetDegree() == 4:
                        # Verify both rings are different
                        if ring1 != ring2:
                            # Additional check: both oxygens should be in different rings
                            if (o1_idx in ring1 and o2_idx in ring2) or (o1_idx in ring2 and o2_idx in ring1):
                                return True, "Contains spiroketal structure: two rings sharing only one carbon atom bonded to two ring oxygens"
    
    return False, "No valid spiroketal structure found"