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
    # Match mono-, di-, tri-, or tetra-phosphate
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2,OX1])[OX2,OX1]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    # Look for ribose/deoxyribose sugar with 5' phosphate
    # More general pattern that matches both ribose and deoxyribose
    # [CH2X4] for 5' carbon, followed by oxygen and phosphorus
    sugar_phosphate = Chem.MolFromSmarts("[CH2X4]OP(~[OX1,OX2])~[OX1,OX2]")
    if not mol.HasSubstructMatch(sugar_phosphate):
        return False, "No 5'-phosphate group found"
        
    # Check for furanose ring with proper substitution pattern
    # More flexible pattern for ribose/deoxyribose
    sugar_ring = Chem.MolFromSmarts("[CH2X4]OP~[#8][C@H]1O[CH]([CH]([CH]1[OX2H0,OX2H1,CH2])[OX2H0,OX2H1])[#7]")
    if not mol.HasSubstructMatch(sugar_ring):
        return False, "No ribose/deoxyribose ring found"

    # Look for nucleobase patterns - more flexible patterns
    # Purine pattern (adenine, guanine and derivatives)
    purine_pattern = Chem.MolFromSmarts("[#7]1[#6]~[#7][#6]2[#6]1~[#7]~[#6]~[#7]2")
    
    # Pyrimidine pattern (cytosine, uracil, thymine and derivatives)
    pyrimidine_pattern = Chem.MolFromSmarts("[#7]1[#6]~[#6]~[#7]~[#6](~[O,S,N])~[#7]1")
    
    has_purine = mol.HasSubstructMatch(purine_pattern)
    has_pyrimidine = mol.HasSubstructMatch(pyrimidine_pattern)
    
    if not (has_purine or has_pyrimidine):
        return False, "No purine or pyrimidine base found"

    # Count phosphate groups
    p_o_p_pattern = Chem.MolFromSmarts("[PX4]~[OX2]~[PX4]")
    p_o_p_count = len(mol.GetSubstructMatches(p_o_p_pattern))
    
    if p_o_p_count > 3:
        return False, "Too many phosphate groups"
        
    phosphate_types = ["mono", "di", "tri", "tetra"]
    phosphate_type = phosphate_types[min(p_o_p_count, 3)]
    
    # Additional check for connection between sugar and base
    sugar_base_connect = Chem.MolFromSmarts("[#7]~[#6,#7]1~[#6,#7]~[#6]~[#6,#7]~1~[CH]1O[CH]([CH2]OP)[CH][CH]1")
    if not mol.HasSubstructMatch(sugar_base_connect):
        return False, "Base not properly connected to sugar"
        
    base_type = "purine" if has_purine else "pyrimidine"
    return True, f"Found {phosphate_type}-phosphorylated nucleoside with {base_type} base"