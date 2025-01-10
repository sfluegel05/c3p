"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: CHEBI:36976 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    A nucleotide must have:
    1. A nucleobase (purine or pyrimidine)
    2. A sugar moiety (typically ribose or deoxyribose)
    3. At least one phosphate group
    4. Correct connectivity between components

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate groups
    phosphate_patterns = [
        "[PX4](=[OX1])([OX2H,OX1-,OX2])([OX2H,OX1-,OX2])[OX2H,OX1-,OX2]",  # Regular phosphate
        "[PX4]([OX2H,OX1-,OX2])([OX2H,OX1-,OX2])([OX2H,OX1-,OX2])[OX2]",   # Phosphate ester
        "[PX4]1([OX2H,OX1-,OX2])(=[OX1])[OX2][CH2][OX2]1",                  # Cyclic phosphate
        "[PX4](=[OX1])([OX2])([OX2H,OX1-,OX2])OP",                          # Di/tri-phosphate start
    ]
    
    has_phosphate = False
    for pattern in phosphate_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_phosphate = True
            break
            
    if not has_phosphate:
        return False, "No phosphate group found"

    # Check for sugar (furanose) ring with more flexible patterns
    sugar_patterns = [
        # Basic furanose patterns (both ribose and deoxyribose)
        "[CH2,CH][OH,O,N]C1OC(C[OH,O,N,P])CC1",
        "[CH2,CH][OH,O,N]C1OC(C[OH,O,N,P])C(O)C1",
        "[CH2,CH][OH,O,N]C1OC(C[OH,O,N,P])C([OH,O])C1",
        # Cyclic patterns
        "C1OC2COP(=O)(O)OC2C1",
        # More general patterns for modified sugars
        "C1OC(CO[P,C])C([OH,O,N,F])C1",
        "C1OC(CO[P,C])CC1",
        # Patterns for various modifications
        "[CH2]1[CH]([OH,O,F,N,S])[CH]([OH,O])[CH]([OH,O])O1",
        # Pattern for N-glycosidic bond connection
        "[$(C1OC(CO)CC1),$(C1OC(CO)C(O)C1)]",
    ]
    
    has_sugar = False
    for pattern in sugar_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_sugar = True
            break
            
    if not has_sugar:
        return False, "No sugar moiety found"

    # Check for nucleobase patterns
    base_patterns = [
        # Purine patterns (adenine, guanine, etc.)
        "c12ncnc([NH2,O])c1nc[nH]2",
        "c12[nH]cnc1c(=O)[nH]c(=O)n2",
        # Pyrimidine patterns (cytosine, uracil, thymine)
        "c1c[nH]c(=O)[nH]c1=O",
        "c1c[nH]c(=O)nc1N",
        # Modified base patterns
        "c1nc([NH2,O])nc2[nH]cnc12",
        "c1nc(N)c2ncn([*])c2n1",
        # General patterns
        "[$(c1nc2c([nH]1)nc[nH]c2=O),$(c1[nH]c(=O)[nH]c(=O)c1)]",
        "[$(c1ncnc2[nH]cnc12),$(c1cncnc1N)]"
    ]
    
    has_base = False
    for pattern in base_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_base = True
            break
            
    if not has_base:
        return False, "No nucleobase found"

    # Check connectivity between components
    connection_patterns = [
        # Sugar-phosphate linkages
        "CO[PX4]",
        "COP(=O)",
        # N-glycosidic bonds
        "c1[nX3]c[nX3]c1[CH]1O[CH]",
        "[nX3]1c[nX3]cc1[CH]1O[CH]",
        # Cyclic nucleotide patterns
        "C1OC2COP(=O)(O)OC2C1",
    ]
    
    has_connections = False
    for pattern in connection_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_connections = True
            break
            
    if not has_connections:
        return False, "Components not properly connected"

    return True, "Contains nucleobase, sugar, and phosphate with correct connectivity"