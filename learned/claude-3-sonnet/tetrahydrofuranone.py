"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
"""
Classifies: tetrahydrofuranone - any oxolane having an oxo- substituent at any position 
on the tetrahydrofuran ring
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone contains a tetrahydrofuran ring with a ketone group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for basic tetrahydrofuran ring with ketone
    # Define patterns for different types of tetrahydrofuranones
    patterns = [
        # Lactone pattern (ketone as part of the ring)
        "[O;R1]1[CH2,CH1][CH2,CH1][CH2,CH1]C1(=O)",
        
        # Ketone at position 2
        "[O;R1]1[C;R1](=O)[CH2,CH1][CH2,CH1][CH2,CH1]1",
        
        # Ketone at position 3
        "[O;R1]1[CH2,CH1][C;R1](=O)[CH2,CH1][CH2,CH1]1",
        
        # Ketone at position 4
        "[O;R1]1[CH2,CH1][CH2,CH1][C;R1](=O)[CH2,CH1]1",
        
        # Exocyclic ketone patterns
        "[O;R1]1[CH2,CH1][CH2,CH1][CH2,CH1][C;R1]1C(=O)",
        "[O;R1]1[CH2,CH1][CH2,CH1][C;R1]([CH2,CH1]1)C(=O)",
        "[O;R1]1[CH2,CH1][C;R1]([CH2,CH1][CH2,CH1]1)C(=O)",
        "[O;R1]1[C;R1]([CH2,CH1][CH2,CH1][CH2,CH1]1)C(=O)"
    ]

    # Check for any of the tetrahydrofuranone patterns
    found_pattern = False
    matched_atoms = set()
    
    for pattern in patterns:
        pattern_mol = Chem.MolFromSmarts(pattern)
        matches = mol.GetSubstructMatches(pattern_mol)
        if matches:
            found_pattern = True
            for match in matches:
                matched_atoms.update(match)
            break

    if not found_pattern:
        return False, "No tetrahydrofuranone core structure found"

    # Check for aromatic character in the matched atoms
    for atom_idx in matched_atoms:
        if mol.GetAtomWithIdx(atom_idx).GetIsAromatic():
            return False, "Core ring structure must be saturated"

    # Count carbonyls
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    
    # Verify that at least one carbonyl is associated with the tetrahydrofuran ring
    carbonyl_connected = False
    for match in carbonyl_matches:
        if any(idx in matched_atoms for idx in match):
            carbonyl_connected = True
            break
            
    if not carbonyl_connected:
        return False, "No ketone group properly connected to tetrahydrofuran ring"

    # Additional structural checks
    ring_info = mol.GetRingInfo()
    ring_systems = ring_info.AtomRings()
    
    # Find the THF ring atoms
    thf_ring = None
    for ring in ring_systems:
        if all(idx in matched_atoms for idx in ring):
            thf_ring = ring
            break
    
    if thf_ring is None:
        return False, "Could not identify clear THF ring"

    # Check if the THF ring is part of a valid structure
    for atom_idx in thf_ring:
        atom = mol.GetAtomWithIdx(atom_idx)
        # Allow normal substituents but prevent certain complex cases
        if atom.GetDegree() > 4:  # Too many connections
            return False, "Invalid substitution pattern on THF ring"

    return True, "Contains tetrahydrofuran ring with ketone group"