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
    A nucleotide must have a nucleoside base, sugar moiety, and at least one phosphate group.

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

    # Check for phosphate groups - expanded patterns
    phosphate_patterns = [
        # Basic phosphate patterns
        Chem.MolFromSmarts("[PX4](=[OX1])([OX2H,OX1-,OX2])([OX2H,OX1-,OX2])[OX2H,OX1-,OX2]"),  # Regular phosphate
        Chem.MolFromSmarts("[PX4]([OX2H,OX1-,OX2])([OX2H,OX1-,OX2])([OX2H,OX1-,OX2])[OX2]"),   # Phosphate ester
        # Di and triphosphates
        Chem.MolFromSmarts("[PX4](=[OX1])([OX2])([OX2H,OX1-,OX2])OP(=[OX1])([OX2H,OX1-,OX2])[OX2H,OX1-,OX2]"),
        Chem.MolFromSmarts("[PX4](=[OX1])([OX2])([OX2H,OX1-,OX2])OP(=[OX1])([OX2])OP(=[OX1])([OX2H,OX1-,OX2])[OX2H,OX1-,OX2]"),
        # Cyclic phosphates
        Chem.MolFromSmarts("[PX4]1([OX2H,OX1-,OX2])(=[OX1])[OX2][CH2][OX2]1"),
        # Mixed phosphate esters
        Chem.MolFromSmarts("[PX4](=[OX1])([OX2C])([OX2H,OX1-,OX2])[OX2H,OX1-,OX2]")
    ]
    
    has_phosphate = any(mol.HasSubstructMatch(pattern) for pattern in phosphate_patterns)
    if not has_phosphate:
        return False, "No phosphate group found"

    # Check for sugar (furanose) ring - expanded patterns
    sugar_patterns = [
        # Ribose patterns
        Chem.MolFromSmarts("[CH2]1[CH]([OH,O])[CH]([OH,O])[CH]([OH,O])O1"),  # Regular ribose
        Chem.MolFromSmarts("[CH2]1[CH]([OH,O])[CH]([OH,O])[CH]O1"),          # Deoxyribose
        # Modified sugar patterns
        Chem.MolFromSmarts("[CH2]1[C@H,C@@H]([OH,O,N])[C@H,C@@H]([OH,O])[C@H,C@@H]([OH,O])O1"),
        Chem.MolFromSmarts("[CH2]1[CH]([OH,O])[CH]([OH,O])[CH]([*])O1"),    # Any substitution
        # Cyclic patterns
        Chem.MolFromSmarts("[CH2]1O[CH]([CH])([CH])[CH]1"),
        # 2'-modified patterns
        Chem.MolFromSmarts("[CH2]1[CH]([OH,O,F,N,S])[CH]([OH,O])[CH]([OH,O])O1")
    ]
    
    has_sugar = any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns)
    if not has_sugar:
        return False, "No sugar (ribose/deoxyribose) moiety found"

    # Check for nucleobase patterns - expanded patterns
    base_patterns = [
        # Purine patterns
        Chem.MolFromSmarts("c12ncnc([NH2,O])c1ncn2"),     # Adenine/Guanine core
        Chem.MolFromSmarts("c12[nH]cnc1c(=O)[nH]c(=O)n2"), # Xanthine core
        # Pyrimidine patterns
        Chem.MolFromSmarts("c1cn([*])c(=O)[nH]c1=O"),     # Uracil/Thymine
        Chem.MolFromSmarts("c1cn([*])c(=O)nc1N"),         # Cytosine
        # Modified base patterns
        Chem.MolFromSmarts("c1nc2c([nH]1)nc[nH]c2=O"),    # Modified purine
        Chem.MolFromSmarts("c1[nX3]c(=O)[nX3]c(=O)c1[*]"), # Modified pyrimidine
        # Additional patterns for rare bases
        Chem.MolFromSmarts("c1nc([NH2,O])nc2c1nc[nH]2"),
        Chem.MolFromSmarts("c1nc(N)c2ncn([*])c2n1"),
        # Methylated/modified bases
        Chem.MolFromSmarts("c1nc(N)c2c1[nH]c[nH]2"),
        Chem.MolFromSmarts("[c,n]1[c,n][c,n][c,n]2[c,n]1[c,n][c,n][nH]2")
    ]
    
    has_base = any(mol.HasSubstructMatch(pattern) for pattern in base_patterns)
    if not has_base:
        return False, "No nucleobase found"

    # Check connectivity between components
    connection_patterns = [
        # Sugar-phosphate linkages
        Chem.MolFromSmarts("[CH2]OP(=[OX1])"),
        Chem.MolFromSmarts("[CH2]OP([OX2H,OX1-,OX2])"),
        # N-glycosidic bonds
        Chem.MolFromSmarts("[NX3]1[CH]([CH2]O)[CH]([OH,O])[CH]([OH,O])[CH]1"),
        Chem.MolFromSmarts("[nX3]1[cX3][nX3][cX3][cX3]1[CH]1O[CH][CH][CH]1"),
        # Cyclic connections
        Chem.MolFromSmarts("[NX3]1[CH]2O[CH][CH][CH]2[CH]([CH2]OP)[CH]1"),
        # Alternative connections
        Chem.MolFromSmarts("[CH2]1[CH]([OH,O])([CH])[CH]([NX3])[O,S]1")
    ]
    
    has_connections = any(mol.HasSubstructMatch(pattern) for pattern in connection_patterns)
    if not has_connections:
        return False, "Components not properly connected"

    return True, "Contains nucleobase, sugar, and phosphate with correct connectivity"