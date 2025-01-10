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
    # More permissive pattern that allows for some variation in bond types
    steroid_core = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1")
    if steroid_core is None:
        return False, "Invalid steroid core SMARTS pattern"
    
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Pattern for 3-oxo group in the first ring
    # This pattern specifically looks for a ketone at position 3 of the A ring
    oxo_pattern = Chem.MolFromSmarts("[#6]~1~[#6]~[#6](=[O])~[#6]~[#6]~[#6]~1~[#6]~2")
    if oxo_pattern is None:
        return False, "Invalid oxo pattern SMARTS"
        
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No ketone group at position 3"

    # Additional validation checks
    
    # Count carbons (steroids typically have 17+ carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 17:
        return False, "Too few carbons for a steroid structure"

    # Check ring count
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings"

    # Verify ring sizes - should have three 6-membered rings and one 5-membered ring
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    six_membered_rings = sum(1 for size in ring_sizes if size == 6)
    five_membered_rings = sum(1 for size in ring_sizes if size == 5)
    
    if six_membered_rings < 3 or five_membered_rings < 1:
        return False, "Incorrect ring sizes for steroid structure"

    # Count oxygen atoms (should have at least one for the ketone)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found"

    # Additional check for ketone specifically
    ketone_pattern = Chem.MolFromSmarts("[#6]-C(=O)-[#6]")
    if ketone_pattern is None:
        return False, "Invalid ketone SMARTS pattern"
    
    ketone_matches = len(mol.GetSubstructMatches(ketone_pattern))
    if ketone_matches < 1:
        return False, "No ketone groups found"

    return True, "Contains steroid core with ketone group at position 3"