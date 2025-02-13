"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
"""
Classifies: CHEBI:16813 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 1-phosphatidyl-1D-myo-inositol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible patterns for myo-inositol ring with correct stereochemistry
    myo_inositol_patterns = [
        # Basic cyclohexane ring with 6 OH groups
        "C1C(O)C(O)C(O)C(O)C(O)C1O",
        # Specific stereochemistry patterns
        "[C@H]1([OH])[C@H]([OH])[C@@H]([OH])[C@H]([OH])[C@@H]([OH])[C@H]1[OH]",
        "[C@@H]1([OH])[C@@H]([OH])[C@H]([OH])[C@@H]([OH])[C@H]([OH])[C@@H]1[OH]",
        # Alternative representations
        "[C@H]1([OH])C([OH])C([OH])C([OH])C([OH])C1[OH]",
        "C1([OH])[C@H]([OH])[C@@H]([OH])[C@H]([OH])[C@@H]([OH])[C@H]1[OH]"
    ]
    
    found_inositol = False
    for pattern in myo_inositol_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt and mol.HasSubstructMatch(patt):
            found_inositol = True
            break
            
    if not found_inositol:
        return False, "No myo-inositol ring found"

    # More flexible phosphate group pattern
    phosphate_patterns = [
        "[O,OH]-P(=O)([O,OH])-O",
        "OP(=O)(O)O",
        "P(=O)(O)(O)O"
    ]
    
    found_phosphate = False
    for pattern in phosphate_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt and mol.HasSubstructMatch(patt):
            found_phosphate = True
            break
            
    if not found_phosphate:
        return False, "No phosphate group found"

    # Check for glycerol backbone - more flexible pattern
    glycerol_patterns = [
        "CC(CO)(CO)",
        "[CH2][CH][CH2]",
        "C[C@H](CO)CO",
        "C[C@@H](CO)CO"
    ]
    
    found_glycerol = False
    for pattern in glycerol_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt and mol.HasSubstructMatch(patt):
            found_glycerol = True
            break
            
    if not found_glycerol:
        return False, "No glycerol backbone found"

    # Check for two ester groups
    ester_pattern = Chem.MolFromSmarts("[#6]-C(=O)-O-[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Check for fatty acid chains (more flexible pattern)
    fatty_acid_pattern = Chem.MolFromSmarts("[#6]~[#6]~[#6]~[#6]~[#6]~[#6]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, "Missing fatty acid chains"

    # Verify overall composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if p_count != 1:
        return False, "Must have exactly one phosphorus atom"
    if o_count < 12:  # 6 from inositol, 2 from esters, 4 from phosphate
        return False, "Insufficient oxygen atoms"
    if c_count < 20:  # minimum for shortest possible chains
        return False, "Carbon chains too short"

    return True, "Contains myo-inositol ring with phosphate group and two fatty acid chains"