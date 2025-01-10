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
    
    # Look for spiro carbon pattern: carbon with two oxygens and two carbons
    # The SMARTS pattern looks for:
    # - A carbon atom [#6] that is
    # - Connected to exactly 4 atoms [X4]
    # - Connected to exactly 2 oxygens <;!@[#8]>{2} (non-aromatic bonds)
    # - Connected to exactly 2 carbons <;!@[#6]>{2} (non-aromatic bonds)
    spiro_pattern = Chem.MolFromSmarts('[#6X4](<;!@[#8]>2)(<;!@[#6]>2)')
    if not mol.HasSubstructMatch(spiro_pattern):
        return False, "No carbon atom with correct connectivity pattern found"

    # Find all spiro carbons
    spiro_matches = mol.GetSubstructMatches(spiro_pattern)
    
    # For each potential spiro carbon
    for match in spiro_matches:
        spiro_carbon = mol.GetAtomWithIdx(match[0])
        
        # Get all neighboring atoms
        neighbors = spiro_carbon.GetNeighbors()
        oxygen_neighbors = [n for n in neighbors if n.GetAtomicNum() == 8]
        carbon_neighbors = [n for n in neighbors if n.GetAtomicNum() == 6]
        
        if len(oxygen_neighbors) != 2 or len(carbon_neighbors) != 2:
            continue
            
        # Check if the oxygens are part of rings
        ring_info = mol.GetRingInfo()
        
        # For each oxygen, follow through its neighbors to find rings
        rings_through_o = []
        for o in oxygen_neighbors:
            for n in o.GetNeighbors():
                if n.GetIdx() != spiro_carbon.GetIdx():
                    rings = ring_info.AtomRings()
                    for ring in rings:
                        if o.GetIdx() in ring and n.GetIdx() in ring:
                            rings_through_o.append(set(ring))
        
        # Need exactly two different rings through the oxygens
        if len(rings_through_o) == 2 and rings_through_o[0] != rings_through_o[1]:
            # Check if the rings share only the spiro carbon
            common_atoms = rings_through_o[0].intersection(rings_through_o[1])
            if len(common_atoms) == 1 and spiro_carbon.GetIdx() in common_atoms:
                return True, "Contains spiroketal structure with two rings sharing only one carbon atom bonded to two oxygens"
    
    return False, "No valid spiroketal structure found"