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
    # More general phosphate pattern that matches mono-, di-, tri- phosphates
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2,OX1])[OX2,OX1]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    # Look for 5'-phosphate connection - CH2-O-P
    sugar_phosphate = Chem.MolFromSmarts("[CH2X4]OP(~[OX1,OX2])~[OX1,OX2]")
    if not mol.HasSubstructMatch(sugar_phosphate):
        return False, "No 5'-phosphate group found"
        
    # Check for furanose ring with more flexible patterns
    sugar_ring_patterns = [
        # Basic furanose ring without stereochemistry
        Chem.MolFromSmarts("[CH2][C]1O[C][C][C]1"),
        # Ribose-like pattern (more general)
        Chem.MolFromSmarts("[CH2]([OH,O])[C]1O[C]([#6,#7])[C][C]1[OH,O]"),
        # Deoxyribose-like pattern (more general)
        Chem.MolFromSmarts("[CH2]([OH,O])[C]1O[C]([#6,#7])[C][C]1"),
        # Alternative sugar pattern
        Chem.MolFromSmarts("[CH2]OP~[OX1,OX2][C]1O[C][C][C]1")
    ]
    
    found_sugar = False
    for pattern in sugar_ring_patterns:
        if mol.HasSubstructMatch(pattern):
            found_sugar = True
            break
    
    if not found_sugar:
        return False, "No furanose ring found"

    # Nucleobase patterns (more general)
    purine_patterns = [
        # Basic purine scaffold (more flexible)
        Chem.MolFromSmarts("[#7,#6]1[#7,#6][#7,#6][#6]2[#7,#6]1[#7,#6][#7,#6][#7,#6]2"),
        # Modified purine
        Chem.MolFromSmarts("[#7,#6]1[#7,#6][#7,#6][#6]2[#7,#6]1[#7,#6][#6]([#7,O,S])[#7,#6]2")
    ]
    
    pyrimidine_patterns = [
        # Basic pyrimidine (more flexible)
        Chem.MolFromSmarts("[#7,#6]1[#7,#6][#7,#6][#6]([#7,O,S])[#7,#6][#7,#6]1"),
        # Alternative pyrimidine
        Chem.MolFromSmarts("[#7]1[#6][#6][#7][#6]([O,S,N])[#7]1")
    ]
    
    has_purine = any(mol.HasSubstructMatch(pattern) for pattern in purine_patterns)
    has_pyrimidine = any(mol.HasSubstructMatch(pattern) for pattern in pyrimidine_patterns)
    
    if not (has_purine or has_pyrimidine):
        return False, "No nucleobase (purine/pyrimidine) found"

    # Count phosphate groups
    p_o_p_pattern = Chem.MolFromSmarts("[PX4]~[OX2]~[PX4]")
    p_o_p_count = len(mol.GetSubstructMatches(p_o_p_pattern))
    phosphate_count = len(mol.GetSubstructMatches(phosphate_pattern))
    
    if p_o_p_count > 3 or phosphate_count > 4:
        return False, "Too many phosphate groups"
        
    # Determine phosphate type
    if p_o_p_count == 0:
        phos_type = "mono"
    elif p_o_p_count == 1:
        phos_type = "di"
    elif p_o_p_count == 2:
        phos_type = "tri"
    else:
        phos_type = "tetra"
    
    base_type = "purine" if has_purine else "pyrimidine"
    return True, f"Found {phos_type}-phosphorylated nucleoside with {base_type} base"