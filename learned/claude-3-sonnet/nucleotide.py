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

    # Check for phosphate group - more general pattern
    phosphate_patterns = [
        Chem.MolFromSmarts("[PX4](=[OX1])([OX2H,OX2-,OX2])[OX2H,OX2-,OX2][OX2H,OX2-,OX2]"),  # Regular phosphate
        Chem.MolFromSmarts("[PX4](=[OX1])([OX2H,OX2-,OX2])[OX2H,OX2-,OX2][OX2]"),  # Phosphoester
        Chem.MolFromSmarts("[PX4]([OX2H,OX2-,OX2])([OX2H,OX2-,OX2])[OX2H,OX2-,OX2][OX2]")  # Cyclic phosphate
    ]
    
    has_phosphate = any(mol.HasSubstructMatch(pattern) for pattern in phosphate_patterns)
    if not has_phosphate:
        return False, "No phosphate group found"

    # Check for sugar (furanose) ring - more general pattern
    sugar_patterns = [
        Chem.MolFromSmarts("[CH2,CH][CH1,CH2][CH1][CH1]O1"),  # Various forms of ribose/deoxyribose
        Chem.MolFromSmarts("[CH2]1O[CH]([CH])[CH]([CH2,CH])[CH]1"),  # Alternative sugar pattern
        Chem.MolFromSmarts("[CH2]1O[CH]([CH,CH2])[CH]([CH,OH])[CH]1")  # More flexible sugar pattern
    ]
    
    has_sugar = any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns)
    if not has_sugar:
        return False, "No sugar (ribose/deoxyribose) moiety found"

    # Check for nucleobase patterns - expanded to catch more variants
    base_patterns = [
        Chem.MolFromSmarts("c12ncnc1[nX3]c[nX3]2"),  # Purine core
        Chem.MolFromSmarts("c1c[nX3]c(=O)[nX3]c1"),  # Pyrimidine core
        Chem.MolFromSmarts("c1nc2c([nX3]1)nc[nX3]2"),  # Alternative purine
        Chem.MolFromSmarts("c1c[nX3]c(=O)nc1[NX3]"),  # Cytosine-like
        Chem.MolFromSmarts("c1[nX3]c(=O)[nX3]c(=O)c1"),  # Uracil/thymine-like
        Chem.MolFromSmarts("c1nc([NX3])nc2c1nc[nX3]2"),  # Modified purine
        Chem.MolFromSmarts("c1nc(O)nc2c1[nX3]c[nX3]2")  # Alternative purine with oxo group
    ]
    
    has_base = any(mol.HasSubstructMatch(pattern) for pattern in base_patterns)
    if not has_base:
        return False, "No nucleobase found"

    # Check for connectivity between components
    # More general patterns for sugar-phosphate and sugar-base connections
    connections = [
        Chem.MolFromSmarts("[CH2,CH]OP(=O)"),  # Sugar-phosphate
        Chem.MolFromSmarts("[CH2,CH]OP([OH,O-])"),  # Alternative sugar-phosphate
        Chem.MolFromSmarts("[NX3]1[CH1][CH1,CH2][CH1][CH1]O1"),  # N-glycosidic bond
        Chem.MolFromSmarts("[NX3][CH]1O[CH][CH][CH]1")  # Alternative N-glycosidic bond
    ]
    
    has_connections = any(mol.HasSubstructMatch(pattern) for pattern in connections)
    if not has_connections:
        return False, "Components not properly connected"

    return True, "Contains nucleobase, sugar, and phosphate with correct connectivity"