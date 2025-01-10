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

    # More flexible steroid core pattern that allows for different bond types and ring variations
    # This pattern matches the four-ring system with more flexibility
    steroid_core = Chem.MolFromSmarts("C1C[C,=C]C2[C,=C]C[C,=C]C3[C,=C][C,=C][C,=C]4[C,=C][C,=C][C,=C]C4(C)[C,=C]3C2(C)C1")
    
    # Alternative core pattern to catch other variations
    steroid_core_alt = Chem.MolFromSmarts("C1CC[C,=C]2[C,=C][C,=C]C3[C,=C][C,=C][C,=C]4[C,=C][C,=C][C,=C]C4(C)C3C2(C)C1")
    
    if not (mol.HasSubstructMatch(steroid_core) or mol.HasSubstructMatch(steroid_core_alt)):
        return False, "No steroid core structure found"

    # Look for ketone at position 3 with more flexible patterns
    # Pattern 1: Standard ketone in ring A
    oxo_pattern1 = Chem.MolFromSmarts("[#6]1[#6][#6](=O)[#6,=C][#6][#6]1")
    # Pattern 2: Alternative ketone arrangement
    oxo_pattern2 = Chem.MolFromSmarts("[#6]1[#6][#6](=O)[#6][#6][#6]1")
    # Pattern 3: Ketone with possible double bonds
    oxo_pattern3 = Chem.MolFromSmarts("[#6]1[#6]C(=O)[#6,=C][#6][#6]1")
    
    if not (mol.HasSubstructMatch(oxo_pattern1) or 
            mol.HasSubstructMatch(oxo_pattern2) or 
            mol.HasSubstructMatch(oxo_pattern3)):
        return False, "No ketone group at position 3"

    # Additional validation checks
    
    # Count carbons (steroids typically have 17+ carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 17:
        return False, "Too few carbons for a steroid structure"

    # Check for proper ring connectivity
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings"

    # Count ketone groups in the molecule
    ketone_pattern = Chem.MolFromSmarts("[#6]-C(=O)-[#6]")
    ketone_matches = len(mol.GetSubstructMatches(ketone_pattern))
    if ketone_matches < 1:
        return False, "No ketone groups found"

    # Check ring sizes (steroids should have mainly 6-membered rings)
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    six_membered_rings = sum(1 for size in ring_sizes if size == 6)
    if six_membered_rings < 2:
        return False, "Not enough 6-membered rings for steroid structure"

    return True, "Contains steroid core with ketone group at position 3"