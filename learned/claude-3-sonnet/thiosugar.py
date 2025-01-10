"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: thiosugar
A carbohydrate derivative where one or more oxygens/hydroxy groups are replaced by sulfur or -SR
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Must contain sulfur
    if not any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms()):
        return False, "No sulfur atoms present"

    # Look for pyranose or furanose ring patterns with more flexible matching
    sugar_patterns = [
        # Pyranose patterns (more flexible)
        "[C]1[O,S][C][C][C][C]1",  # Basic ring
        "[C]1[O][C]([!O])[C][C][C]1", # Ring with substitution at C1
        # Furanose patterns
        "[C]1[O,S][C][C][C]1",
        "[C]1[O][C]([!O])[C][C]1"
    ]
    
    found_sugar = False
    for pattern in sugar_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_sugar = True
            break
            
    if not found_sugar:
        return False, "No sugar-like ring structure found"

    # Check for various sulfur-containing patterns
    s_patterns = [
        "[C]S[C,O,N]",  # S-glycosidic or similar
        "[C]1[O][C]([S])[C,O,N]", # S at anomeric position
        "[C]1[O][C][C]([S])[C]", # S substitution on ring
        "[C]S(=O)", # Sulfoxide
        "[C]S(=O)(=O)", # Sulfone
        "OS(=O)(=O)[O-,OH]" # Sulfonate
    ]
    
    found_s_pattern = False
    for pattern in s_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_s_pattern = True
            break
            
    if not found_s_pattern:
        return False, "No relevant sulfur-containing groups found"

    # Count carbons and heteroatoms in the core structure
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    
    # Basic checks for sugar-like composition
    if c_count < 4:
        return False, "Too few carbons for a sugar structure"
    
    if o_count + s_count < 3:
        return False, "Too few O/S atoms for a sugar structure"

    # Look for potential substituents that are common in thiosugars
    substituent_patterns = [
        "OC(=O)[C]", # Ester
        "S[C](=N)", # Thioglycoside with imine
        "S[C](=O)", # Thioester
        "SC", # Simple thioether
        "OS(=O)(=O)[O-,OH]" # Sulfonate
    ]
    
    has_valid_substituent = False
    for pattern in substituent_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_valid_substituent = True
            break

    # Final check combining all evidence
    if found_sugar and found_s_pattern and (has_valid_substituent or (o_count + s_count >= 4)):
        return True, "Contains sugar ring with sulfur-containing modification"
    
    return False, "Structure does not match thiosugar patterns"