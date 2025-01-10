"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
"""
Classifies: nucleoside phosphate
A nucleobase-containing molecular entity that is a nucleoside in which one or more 
of the sugar hydroxy groups has been converted into a mono- or poly-phosphate.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a nucleoside phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for phosphate groups - multiple patterns to catch variations
    phosphate_patterns = [
        Chem.MolFromSmarts("[OX2,OX1-]P(=O)([OX2,OX1-])[OX2,OX1-]"),  # Standard phosphate
        Chem.MolFromSmarts("[OH]P(=O)([OH])[OH]"),  # Fully protonated
        Chem.MolFromSmarts("[O-]P(=O)([O-])[O-]"),  # Fully deprotonated
        Chem.MolFromSmarts("OP(=O)(O)OP(=O)(O)O"),  # Diphosphate
        Chem.MolFromSmarts("OP(=O)(O)OP(=O)(O)OP(=O)(O)O"),  # Triphosphate
    ]
    
    has_phosphate = any(mol.HasSubstructMatch(pattern) for pattern in phosphate_patterns if pattern is not None)
    if not has_phosphate:
        return False, "No phosphate group found"
    
    # Check for nucleobase patterns - expanded set
    purine_patterns = [
        Chem.MolFromSmarts("c1ncnc2[nX3]cnc12"),  # Adenine core
        Chem.MolFromSmarts("c1nc(N)nc2[nX3]cnc12"),  # Adenine
        Chem.MolFromSmarts("c1nc(=O)[nH]c2[nX3]cnc12"),  # Guanine core
        Chem.MolFromSmarts("c1nc(N)c2nc[nH]c2n1"),  # Alternative purine
        Chem.MolFromSmarts("[nX3]1cnc2c(ncnc2[nX3]1)N"),  # Modified purine
    ]
    
    pyrimidine_patterns = [
        Chem.MolFromSmarts("[nX3]1c([nX3H])cc(=O)[nX3]c1=O"),  # Basic pyrimidine
        Chem.MolFromSmarts("O=c1cc[nH]c(=O)[nX3]1"),  # Uracil
        Chem.MolFromSmarts("O=c1cc[nH]c(=O)n1[CH3]"),  # Thymine
        Chem.MolFromSmarts("Nc1cc[nX3]c(=O)[nX3]1"),  # Cytosine
        Chem.MolFromSmarts("[nX3]1ccc(=O)[nX3]c1=O"),  # Alternative pyrimidine
    ]
    
    has_purine = any(mol.HasSubstructMatch(pattern) for pattern in purine_patterns if pattern is not None)
    has_pyrimidine = any(mol.HasSubstructMatch(pattern) for pattern in pyrimidine_patterns if pattern is not None)
    
    if not (has_purine or has_pyrimidine):
        return False, "No nucleobase (purine or pyrimidine) found"
    
    # Check for sugar (ribose/deoxyribose) patterns - more flexible
    sugar_patterns = [
        Chem.MolFromSmarts("[CH2X4]-[CH1X4]-[CH1X4]-[CH1X4]-[OX2]"),  # Basic sugar
        Chem.MolFromSmarts("[CH2X4]-[CX4]-[CX4]-[CH1X4]-[OX2]"),  # Deoxyribose
        Chem.MolFromSmarts("[CH2X4]-1-[CH1X4]-[CH1X4]-[CH1X4]-O1"),  # Cyclic sugar
    ]
    
    has_sugar = any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns if pattern is not None)
    if not has_sugar:
        return False, "No sugar moiety found"
    
    # Check connectivity between sugar and phosphate - more permissive
    sugar_phosphate_patterns = [
        Chem.MolFromSmarts("[CH2X4]OP(=O)([O,OH,O-])[O,OH,O-]"),  # 5' phosphate
        Chem.MolFromSmarts("[CH1X4]OP(=O)([O,OH,O-])[O,OH,O-]"),  # 2' or 3' phosphate
        Chem.MolFromSmarts("[CH2X4]OP([O,OH,O-])([O,OH,O-])OP"),  # Di/tri-phosphate
    ]
    
    has_connection = any(mol.HasSubstructMatch(pattern) for pattern in sugar_phosphate_patterns if pattern is not None)
    if not has_connection:
        return False, "Sugar and phosphate not properly connected"
    
    # Count phosphates for reporting
    phosphate_count = sum(len(mol.GetSubstructMatches(pattern)) 
                         for pattern in phosphate_patterns if pattern is not None)
    base_type = "purine" if has_purine else "pyrimidine"
    
    return True, f"Contains {base_type} nucleobase, sugar moiety, and {phosphate_count} phosphate group(s)"