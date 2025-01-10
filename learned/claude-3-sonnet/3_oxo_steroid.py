"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: CHEBI:35341 3-oxo steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid has a steroid core structure with a ketone group at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic steroid core pattern - four fused rings (three 6-membered, one 5-membered)
    # More flexible pattern that allows for different bond types and substitutions
    steroid_core = Chem.MolFromSmarts(
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6,#8]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1"
    )
    
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Pattern for 3-oxo group
    # This looks specifically for a ketone at position 3 in the A ring
    # The pattern ensures it's connected to the steroid core properly
    oxo_pattern = Chem.MolFromSmarts(
        "[#6]~1~[#6]~[#6](=O)~[#6]~[#6]~[#6]~1~[#6]~2~[#6]~[#6]~[#6]"
    )
    
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No ketone group at position 3"

    # Count carbons (steroids typically have 17+ carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 17:
        return False, "Too few carbons for a steroid structure"

    # Check ring count and sizes
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    
    # Should have at least 4 rings
    if len(ring_sizes) < 4:
        return False, "Insufficient number of rings"

    # Count 5 and 6 membered rings
    six_membered = sum(1 for size in ring_sizes if size == 6)
    five_membered = sum(1 for size in ring_sizes if size == 5)
    
    if six_membered < 2 or five_membered < 1:
        return False, "Incorrect ring composition for steroid structure"

    # Verify presence of ketone group
    ketone_pattern = Chem.MolFromSmarts("[#6]-C(=O)-[#6]")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group found"

    # More specific check for 3-oxo position using a detailed pattern
    # This pattern ensures the ketone is at position 3 relative to the ring fusion
    specific_3_oxo = Chem.MolFromSmarts(
        "[#6]~1~[#6]~[#6](=O)~[#6]~[#6]~[#6]~1~[#6]~2~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]"
    )
    
    if not mol.HasSubstructMatch(specific_3_oxo):
        return False, "Ketone group not at position 3"

    return True, "Contains steroid core with ketone group at position 3"