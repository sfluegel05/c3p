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

    # First check for the basic steroid skeleton with more flexible patterns
    # This pattern allows for variations in bond types and substitution
    steroid_pattern = Chem.MolFromSmarts(
        "[C,c]1~[C,c]~[C,c]~[C,c]2~[C,c]~[C,c]~[C,c]3~[C,c]~[C,c]~[C,c]4~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]4~[C,c]~[C,c]3~[C,c]~[C,c]2~[C,c]1"
    )
    
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid core structure found"

    # Check for the specific 3-oxo-Delta(4) pattern in ring A
    # This pattern looks for:
    # - A ketone group at position 3
    # - A double bond between positions 4 and 5
    oxo_delta4_pattern = Chem.MolFromSmarts(
        "[C,c]1~[C,c]~C(=O)~C=C~[C,c]2~[C,c]~[C,c]~[C,c]1~[C,c]2"
    )
    
    if not mol.HasSubstructMatch(oxo_delta4_pattern):
        return False, "No 3-oxo-Delta(4) pattern found"

    # Additional validation
    # Count carbons (steroids typically have 19-35 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (19 <= carbon_count <= 35):
        return False, f"Invalid carbon count ({carbon_count}) for steroid structure"

    # Count rings (steroids must have at least 4 rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Check for ketone group specifically
    ketone_pattern = Chem.MolFromSmarts("[#6]-C(=O)-[#6]")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group found"

    # Check for conjugated double bond
    conjugated_db_pattern = Chem.MolFromSmarts("C(=O)-C=C")
    if not mol.HasSubstructMatch(conjugated_db_pattern):
        return False, "No conjugated double bond found"

    return True, "Contains steroid core with 3-oxo group conjugated to Delta-4 double bond"