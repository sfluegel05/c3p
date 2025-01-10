"""
Classifies: CHEBI:18133 hexose
"""
"""
Classifies: CHEBI:24895 hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Pyranose patterns (6-membered ring)
    pyranose_patterns = [
        # Basic pyranose with various substitutions
        "[C]-1-[C]-[C]-[C]-[C](-[O]-1)",
        # Alternative forms with different oxidation states
        "O1[C][C][C][C][C]1",
        "O1[C][C][C][C](O)[C]1"
    ]
    
    # Furanose patterns (5-membered ring)
    furanose_patterns = [
        # Basic furanose with various substitutions
        "[C]-1-[C]-[C]-[C](-[O]-1)-[C]",
        # Alternative forms
        "O1[C][C][C][C]1[C]",
        "O1[C][C][C][C](CO)1"
    ]
    
    # Open chain patterns
    open_chain_patterns = [
        # Aldehexoses
        "[CH](=O)[C][C][C][C][C]",
        # Ketohexoses
        "[C][C](=O)[C][C][C][C]",
        # Alternative forms
        "[CH](=O)[C][C][C][C]CO",
        "[C][C](=O)[C][C][C]CO"
    ]

    # Check for any hexose core structure
    found_core = False
    for pattern in pyranose_patterns + furanose_patterns + open_chain_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt and mol.HasSubstructMatch(patt):
            found_core = True
            break
    
    if not found_core:
        return False, "No hexose core structure found"

    # Count carbons (should be around 6, but allow for modifications)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 6:
        return False, f"Too few carbons for hexose ({carbon_count})"
    if carbon_count > 12:  # Allow for some substitutions
        return False, f"Too many carbons for simple hexose ({carbon_count})"

    # Check for oxygen content
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 5:  # Most hexoses have at least 5 oxygens (ring + hydroxyls)
        return False, f"Too few oxygens for hexose ({oxygen_count})"

    # Check for characteristic hydroxyl/carbonyl patterns
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    
    has_hydroxyls = mol.HasSubstructMatch(hydroxyl_pattern)
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)
    
    if not (has_hydroxyls or has_carbonyl):
        return False, "Missing characteristic hydroxyl/carbonyl groups"

    # Check it's not a disaccharide
    glycosidic_pattern = Chem.MolFromSmarts("[OX2]([CH1][OX2])[CH1]")
    if mol.HasSubstructMatch(glycosidic_pattern):
        matches = len(mol.GetSubstructMatches(glycosidic_pattern))
        if matches > 1:
            return False, "Appears to be a polysaccharide"

    # Additional validation for carbon chain connectivity
    chain_patterns = [
        "[C]-[C]-[C]-[C]-[C]-[C]",
        "[C]1[C][C][C][C][C]1",
        "[C]1[C][C][C][C]1[C]"
    ]
    
    chain_found = False
    for pattern in chain_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt and mol.HasSubstructMatch(patt):
            chain_found = True
            break
            
    if not chain_found:
        return False, "Missing connected six-carbon chain"

    return True, "Contains hexose structure with appropriate substitution pattern"