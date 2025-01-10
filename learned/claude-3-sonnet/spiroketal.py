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
    # [#6X4] is sp3 carbon
    # [#8X2;R] is ring oxygen with 2 connections
    # Pattern looks for sp3 carbon with two single-bonded ring oxygens
    spiro_pattern = Chem.MolFromSmarts('[#6X4](-[#8X2;R])-[#8X2;R]')
    
    if not mol.HasSubstructMatch(spiro_pattern):
        return False, "No sp3 carbon atom with two ring oxygens found"
    
    matches = mol.GetSubstructMatches(spiro_pattern)
    
    # For each potential spiro carbon
    for match in matches:
        spiro_carbon_idx = match[0]
        o1_idx = match[1]
        o2_idx = match[2]
        
        # Get ring information
        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()
        
        # Find smallest rings containing each oxygen
        rings_o1 = []
        rings_o2 = []
        
        for ring in rings:
            if o1_idx in ring:
                rings_o1.append(set(ring))
            if o2_idx in ring:
                rings_o2.append(set(ring))
                
        # If either oxygen isn't in any ring, skip
        if not rings_o1 or not rings_o2:
            continue
            
        # Get smallest rings for each oxygen
        smallest_ring_o1 = min(rings_o1, key=len)
        smallest_ring_o2 = min(rings_o2, key=len)
        
        # Check if rings are different and share only the spiro carbon
        intersection = smallest_ring_o1.intersection(smallest_ring_o2)
        if len(intersection) == 1 and spiro_carbon_idx in intersection:
            # Verify both rings are of reasonable size (3-8 members is typical)
            if 3 <= len(smallest_ring_o1) <= 8 and 3 <= len(smallest_ring_o2) <= 8:
                # Verify the spiro carbon is truly tetrahedral
                spiro_carbon = mol.GetAtomWithIdx(spiro_carbon_idx)
                if spiro_carbon.GetHybridization() == Chem.HybridizationType.SP3:
                    # Verify oxygens are in different rings
                    if (o1_idx in smallest_ring_o1 and o2_idx in smallest_ring_o2) or \
                       (o1_idx in smallest_ring_o2 and o2_idx in smallest_ring_o1):
                        # Verify bonds to oxygens are single bonds
                        bond1 = mol.GetBondBetweenAtoms(spiro_carbon_idx, o1_idx)
                        bond2 = mol.GetBondBetweenAtoms(spiro_carbon_idx, o2_idx)
                        if bond1.GetBondType() == Chem.BondType.SINGLE and \
                           bond2.GetBondType() == Chem.BondType.SINGLE:
                            return True, "Contains spiroketal structure: two rings sharing only one carbon atom bonded to two ring oxygens"
    
    return False, "No valid spiroketal structure found"