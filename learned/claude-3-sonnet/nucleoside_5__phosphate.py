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
    
    # Look for phosphate group (-OP(=O)(O)-)
    phosphate_pattern = Chem.MolFromSmarts("[OX2][PX4](=[OX1])([OX2])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    # Look for ribose or deoxyribose sugar
    # Pattern matches both ribose and deoxyribose with a phosphate at 5' position
    sugar_pattern = Chem.MolFromSmarts("[CH2X4]OP([OX2,OX1])~[#8][C@H]1O[CH]([CH]([CH]1[OX2,CH2])[OX2H0,OX2H1])[#7]")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No ribose/deoxyribose sugar with 5' phosphate found"
    
    # Look for nucleobase patterns
    # Purine pattern (adenine, guanine and derivatives)
    purine_pattern = Chem.MolFromSmarts("c12ncnc([NH2,O])c1nc([NH2,O])[nH]2")
    # Pyrimidine pattern (cytosine, uracil, thymine and derivatives)
    pyrimidine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#7][#6](=[O,S])[#7]1")
    
    has_purine = mol.HasSubstructMatch(purine_pattern)
    has_pyrimidine = mol.HasSubstructMatch(pyrimidine_pattern)
    
    if not (has_purine or has_pyrimidine):
        return False, "No purine or pyrimidine base found"
    
    # Count phosphate groups (look for P-O-P bonds for di/tri/tetra-phosphates)
    p_o_p_pattern = Chem.MolFromSmarts("[PX4]-[OX2]-[PX4]")
    p_o_p_matches = len(mol.GetSubstructMatches(p_o_p_pattern))
    
    phosphate_type = "mono"
    if p_o_p_matches == 1:
        phosphate_type = "di"
    elif p_o_p_matches == 2:
        phosphate_type = "tri"
    elif p_o_p_matches == 3:
        phosphate_type = "tetra"
    
    # Verify connection between sugar and base
    sugar_base_pattern = Chem.MolFromSmarts("[#7]1[#6,#7][#6,#7][#6][#6,#7]1-[CH]1O[CH]([CH2]OP)[CH][CH]1")
    if not mol.HasSubstructMatch(sugar_base_pattern):
        return False, "Base not properly connected to sugar"
    
    base_type = "purine" if has_purine else "pyrimidine"
    return True, f"Found {phosphate_type}-phosphorylated nucleoside with {base_type} base"