"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
"""
Classifies: tetrahydrofuranone - any oxolane having an oxo- substituent at any position 
on the tetrahydrofuran ring
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import GetSymmSSSR

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

    # Look for basic tetrahydrofuran ring
    thf_pattern = Chem.MolFromSmarts("[O;R1]1[CH2,CH1,C][CH2,CH1,C][CH2,CH1,C][CH2,CH1,C]1")
    if not mol.HasSubstructMatch(thf_pattern):
        return False, "No tetrahydrofuran ring found"

    # Check if the ring is aromatic
    aromatic_thf = Chem.MolFromSmarts("[o]1cccc1")
    if mol.HasSubstructMatch(aromatic_thf):
        return False, "Core ring structure must be saturated (tetrahydro-)"

    # Count rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 2:  # Allow for at most 2 rings (THF + possibly one more)
        return False, "Structure too complex - too many rings"

    # Define patterns for different types of tetrahydrofuranones
    patterns = [
        # Simple lactone (ketone as part of the ring)
        Chem.MolFromSmarts("O1CC[CH2,CH1,C]C(=O)1"),
        
        # Ketone at position 2
        Chem.MolFromSmarts("[O;R1]1[C;R1](=O)[CH2,CH1,C][CH2,CH1,C][CH2,CH1,C]1"),
        
        # Ketone at position 3
        Chem.MolFromSmarts("[O;R1]1[CH2,CH1,C][C;R1](=O)[CH2,CH1,C][CH2,CH1,C]1"),
        
        # Ketone at position 4
        Chem.MolFromSmarts("[O;R1]1[CH2,CH1,C][CH2,CH1,C][C;R1](=O)[CH2,CH1,C]1"),
        
        # Ketone substituent patterns
        Chem.MolFromSmarts("[O;R1]1[CH2,CH1,C][CH2,CH1,C][CH2,CH1,C]([CH2,CH1,C]1)C(=O)"),
        Chem.MolFromSmarts("[O;R1]1[CH2,CH1,C][CH2,CH1,C][CH2,CH1,C]([C;R1]1)C(=O)"),
    ]

    # Check for any of the tetrahydrofuranone patterns
    found_pattern = False
    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            found_pattern = True
            break

    if not found_pattern:
        return False, "No ketone group properly positioned on tetrahydrofuran ring"

    # Additional check for ring fusion
    ring_atoms = set()
    for ring in ring_info.AtomRings():
        ring_atoms.update(ring)
    if len(ring_atoms) < sum(len(ring) for ring in ring_info.AtomRings()):
        # If rings share atoms (fused), be more selective
        if ring_info.NumRings() > 1:  # If there's more than one ring and they're fused
            return False, "Complex fused ring system detected"

    return True, "Contains tetrahydrofuran ring with ketone group"