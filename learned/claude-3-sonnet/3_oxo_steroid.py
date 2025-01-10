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

    # Basic steroid core pattern (four fused rings)
    steroid_core = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Look for ketone at position 3 
    # The pattern looks for a ketone group connected to the A ring of the steroid
    oxo_pattern = Chem.MolFromSmarts("[#6]1[#6][#6](=O)[#6][#6][#6]1")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No ketone group at position 3"

    # Additional checks to confirm steroid nature:
    
    # Count carbons (steroids typically have 17+ carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 17:
        return False, "Too few carbons for a steroid structure"

    # Check for proper ring connectivity
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings"

    # Look for common modifications that would invalidate the structure
    # (none identified that would specifically invalidate a 3-oxo steroid)

    return True, "Contains steroid core with ketone group at position 3"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35341',
        'name': '3-oxo steroid',
        'definition': 'Any oxo steroid where an oxo substituent is located at position 3.',
        'parents': ['CHEBI:35350']
    }
}