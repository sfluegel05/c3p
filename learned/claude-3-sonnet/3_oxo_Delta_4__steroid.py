"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
"""
Classifies: CHEBI:35353 3-oxo-Delta(4) steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    These steroids have a ketone at position 3 conjugated with a C=C double bond at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(4) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic steroid core (four fused rings)
    # More flexible pattern that focuses on connectivity rather than specific bond types
    steroid_core = Chem.MolFromSmarts(
        "[#6]1[#6][#6][#6]2[#6][#6][#6]3[#6][#6][#6]4[#6][#6][#6][#6][#6]4[#6][#6]3[#6][#6]2[#6]1"
    )
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for 3-oxo-Delta(4) pattern in A-ring
    # This pattern specifically looks for:
    # - A ketone at position 3
    # - A double bond between positions 4 and 5
    # - Proper connectivity to the B ring
    oxo_delta4_pattern = Chem.MolFromSmarts(
        "[#6]1[#6][#6](=O)[#6]=[#6][#6]2[#6][#6][#6]1[#6]2"
    )
    if not mol.HasSubstructMatch(oxo_delta4_pattern):
        return False, "No 3-oxo-Delta(4) pattern found in A-ring"

    # Additional validation for steroid structure
    # Count carbons (steroids typically have 19-30 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 19 or carbon_count > 35:
        return False, f"Invalid carbon count ({carbon_count}) for steroid structure"

    # Check ring count (steroids must have exactly 4 rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Check for sp3 carbons at ring junctions (characteristic of steroids)
    ring_junction_pattern = Chem.MolFromSmarts("[#6D4]1[#6][#6][#6]2")
    if not mol.HasSubstructMatch(ring_junction_pattern):
        return False, "Missing characteristic ring junction pattern"

    # Verify ketone is at position 3 by checking its environment
    ketone_pattern = Chem.MolFromSmarts(
        "[#6]1[#6][#6](=O)[#6]=[#6][#6]2[#6][#6][#6]1[#6]2"
    )
    if len(mol.GetSubstructMatches(ketone_pattern)) != 1:
        return False, "Ketone not in correct position or multiple ketones found"

    return True, "Contains steroid core with 3-oxo group conjugated to Delta-4 double bond"