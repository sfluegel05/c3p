"""
Classifies: CHEBI:18133 hexose
"""
"""
Classifies: CHEBI:24895 hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import ReactionFromSmarts

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is a six-carbon monosaccharide with either an aldehyde group (aldohexose)
    or a ketone group (ketohexose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define hexose core patterns
    # Pyranose patterns (both alpha and beta)
    pyranose_patterns = [
        "[C]-1-[C]-[C]-[C]-[C](-[O]-1)(-[OH])-[CH2][OH]",  # Basic pyranose
        "[C]-1-[C]-[C]-[C]-[C](-[O]-1)-[CH2][OH]",         # Alternative form
        "[C]-1-[C]-[C]-[C]-[C](-[O]-1)-[CH2]-[O]"          # Modified form
    ]
    
    # Furanose patterns
    furanose_patterns = [
        "[C]-1-[C]-[C]-[C](-[O]-1)-[C](-[OH])-[CH2][OH]",  # Basic furanose
        "[C]-1-[C]-[C]-[C](-[O]-1)-[C]-[CH2][OH]",         # Alternative form
        "[C]-1-[C]-[C]-[C](-[O]-1)-[C]-[CH2]-[O]"          # Modified form
    ]
    
    # Open chain patterns (aldehexoses and ketohexoses)
    open_chain_patterns = [
        "[CH](=O)-[C]-[C]-[C]-[C]-[CH2][OH]",              # Aldehexose
        "[C](=O)-[C]-[C]-[C]-[C]-[CH2][OH]",               # Ketohexose
        "[CH](=O)-[C]-[C]-[C]-[C]-[CH2]-[O]"               # Modified form
    ]

    # Check for hexose core structure
    found_core = False
    for pattern in pyranose_patterns + furanose_patterns + open_chain_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(patt):
            found_core = True
            break
    
    if not found_core:
        return False, "No hexose core structure found"

    # Count carbons in the core structure (should be 6)
    core_carbons = len(mol.GetSubstructMatch(patt))
    if core_carbons < 6:
        return False, f"Core structure contains fewer than 6 carbons ({core_carbons})"

    # Check for characteristic hydroxyl pattern
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    if mol.HasSubstructMatch(hydroxyl_pattern):
        hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
        if hydroxyl_count < 2:  # Allow for modified hexoses
            return False, f"Too few hydroxyl groups ({hydroxyl_count})"
    
    # Check it's not a disaccharide or larger
    glycosidic_pattern = Chem.MolFromSmarts("[OX2]([CH1][OX2])[CH1]")
    if mol.HasSubstructMatch(glycosidic_pattern):
        glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
        if len(glycosidic_matches) > 1:
            return False, "Appears to be a polysaccharide"

    # Additional check for characteristic carbon chain
    chain_pattern = Chem.MolFromSmarts("[C]-[C]-[C]-[C]-[C]-[C]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Missing characteristic six-carbon chain"

    return True, "Contains hexose core structure with appropriate substitution pattern"