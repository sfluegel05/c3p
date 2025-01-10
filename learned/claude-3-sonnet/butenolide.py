"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: CHEBI:35631 butenolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide is a gamma-lactone containing a 2-furanone skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for butenolide core
    # Pattern 1: Basic 2-furanone skeleton with double bond
    # O=C1OCC=C1
    butenolide_pattern1 = Chem.MolFromSmarts("O=C1OCC=C1")
    
    # Pattern 2: Alternative pattern with double bond in different position
    # O=C1OC=CC1
    butenolide_pattern2 = Chem.MolFromSmarts("O=C1OC=CC1")
    
    # Check for basic butenolide core
    has_pattern1 = mol.HasSubstructMatch(butenolide_pattern1)
    has_pattern2 = mol.HasSubstructMatch(butenolide_pattern2)
    
    if not (has_pattern1 or has_pattern2):
        return False, "No 2-furanone core found"

    # Additional checks to exclude false positives
    
    # Count ring size where the pattern was found
    ring_info = mol.GetRingInfo()
    if not any(len(ring) == 5 for ring in ring_info.AtomRings()):
        return False, "No 5-membered ring found"
    
    # Verify lactone group (cyclic ester)
    lactone_pattern = Chem.MolFromSmarts("[O;R1]1[#6][#6](=[O;R0])[O;R1]1")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone group found"

    # Count number of double bonds in the ring
    db_pattern = Chem.MolFromSmarts("C=C")
    db_matches = len(mol.GetSubstructMatches(db_pattern))
    if db_matches == 0:
        return False, "No double bond found in structure"

    # Success - molecule contains butenolide core
    if has_pattern1:
        return True, "Contains 2-furanone core with C3-C4 double bond"
    else:
        return True, "Contains 2-furanone core with C4-C5 double bond"