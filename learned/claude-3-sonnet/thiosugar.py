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

    # Define and validate sugar ring patterns
    pyranose_pattern = "[C,O,S]1[C,O,S][C,O,S][C,O,S][C,O,S][C,O,S]1" # 6-membered ring
    furanose_pattern = "[C,O,S]1[C,O,S][C,O,S][C,O,S][C,O,S]1" # 5-membered ring
    
    ring_patterns = [pyranose_pattern, furanose_pattern]
    found_sugar_ring = False
    
    for pattern in ring_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is not None and mol.HasSubstructMatch(patt):
            found_sugar_ring = True
            break
            
    if not found_sugar_ring:
        return False, "No sugar ring structure found"

    # Check for hydroxyl groups on ring carbons (typical for sugars)
    hydroxyl_pattern = Chem.MolFromSmarts("[C]1[C,O,S][C,O,S][C,O,S][C,O,S][C,O,S]1[OH]")
    if hydroxyl_pattern is None or not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing characteristic sugar hydroxyl groups"

    # Check for various sulfur-containing patterns
    s_patterns = [
        "CS[C,O,N]",  # Thioglycosidic or similar
        "CS(=O)",  # Sulfoxide
        "CS(=O)(=O)", # Sulfone
        "COS(=O)(=O)O", # Sulfonate
        "CS[C@@H]", # Chiral thioglycoside
        "CS[C@H]"  # Chiral thioglycoside
    ]
    
    found_s_pattern = False
    for pattern in s_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is not None and mol.HasSubstructMatch(patt):
            found_s_pattern = True
            break
            
    if not found_s_pattern:
        return False, "No characteristic sulfur-containing groups found"

    # Count key atoms and check ratios
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    
    # Basic composition checks
    if c_count < 4:
        return False, "Too few carbons for a sugar structure"
    
    if o_count < 2:
        return False, "Too few oxygens for a sugar structure"
        
    if o_count + s_count < 4:
        return False, "Too few O/S atoms for a thiosugar"

    # Check for chiral centers (sugars typically have multiple)
    chiral_centers = len(Chem.FindMolChiralCenters(mol))
    if chiral_centers < 2:
        return False, "Too few chiral centers for a sugar derivative"

    # Final check combining all evidence
    if found_sugar_ring and found_s_pattern:
        reason = "Contains sugar ring with "
        if s_count > 1:
            reason += f"{s_count} sulfur-containing modifications"
        else:
            reason += "sulfur-containing modification"
        return True, reason
    
    return False, "Structure does not match thiosugar patterns"