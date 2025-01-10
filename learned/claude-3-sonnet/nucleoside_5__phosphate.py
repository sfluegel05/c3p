"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: CHEBI:24836 nucleoside 5'-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for phosphate group(s)
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2,OX1])[OX2,OX1]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    # Look for 5'-phosphate connection
    # CH2-O-P where CH2 is part of sugar ring
    sugar_phosphate = Chem.MolFromSmarts("[CH2X4]OP(~[OX1,OX2])~[OX1,OX2]")
    if not mol.HasSubstructMatch(sugar_phosphate):
        return False, "No 5'-phosphate group found"
        
    # Check for furanose ring - more flexible pattern
    # Match both ribose and deoxyribose with different possible substituents
    sugar_ring_patterns = [
        # Ribose pattern
        Chem.MolFromSmarts("[CH2X4]O[P]~[#8][C@H]1O[CH]([CH]([CH]1O)[OH])[#7,#6]"),
        # Deoxyribose pattern
        Chem.MolFromSmarts("[CH2X4]O[P]~[#8][C@H]1O[CH]([CH]([CH]1)[H])[#7,#6]")
    ]
    
    found_sugar = False
    for pattern in sugar_ring_patterns:
        if mol.HasSubstructMatch(pattern):
            found_sugar = True
            break
    
    if not found_sugar:
        return False, "No ribose/deoxyribose ring found"

    # Look for nucleobase patterns
    # Purine pattern (adenine, guanine and derivatives)
    purine_patterns = [
        # Basic purine scaffold
        Chem.MolFromSmarts("[#7]1[#6]~[#7][#6]2[#6]1~[#7]~[#6]~[#7]2"),
        # Modified purine scaffold
        Chem.MolFromSmarts("[#7]1[#6]~[#7][#6]2[#6]1~[#7]~[#6](~[O,N])~[#7]2")
    ]
    
    # Pyrimidine pattern (cytosine, uracil, thymine and derivatives)
    pyrimidine_patterns = [
        # Basic pyrimidine scaffold
        Chem.MolFromSmarts("[#7]1[#6]~[#6]~[#7]~[#6](~[O,S,N])~[#7]1"),
        # Modified pyrimidine scaffold
        Chem.MolFromSmarts("[#7]1[#6](~[O,N])[#7][#6]~[#6]~[#7]1")
    ]
    
    has_purine = any(mol.HasSubstructMatch(pattern) for pattern in purine_patterns)
    has_pyrimidine = any(mol.HasSubstructMatch(pattern) for pattern in pyrimidine_patterns)
    
    if not (has_purine or has_pyrimidine):
        return False, "No purine or pyrimidine base found"

    # Count phosphate groups to determine type
    p_o_p_pattern = Chem.MolFromSmarts("[PX4]~[OX2]~[PX4]")
    p_o_p_count = len(mol.GetSubstructMatches(p_o_p_pattern))
    
    if p_o_p_count > 3:
        return False, "Too many phosphate groups"
        
    phosphate_types = ["mono", "di", "tri", "tetra"]
    phosphate_type = phosphate_types[min(p_o_p_count, 3)]
    
    # Check for proper connection between sugar and base
    # More flexible pattern that allows for various modifications
    sugar_base_connects = [
        Chem.MolFromSmarts("[#7]~[#6,#7]1~[#6,#7]~[#6]~[#6,#7]~1~[CH]1O[CH]([CH2]OP)[CH][CH]1"),
        Chem.MolFromSmarts("[#7]~[#6,#7]1~[#6,#7]~[#6]~[#6,#7]~1~[CH]1O[CH]([CH2]OP)[CH][CH]1")
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in sugar_base_connects):
        return False, "Base not properly connected to sugar"
        
    base_type = "purine" if has_purine else "pyrimidine"
    return True, f"Found {phosphate_type}-phosphorylated nucleoside with {base_type} base"